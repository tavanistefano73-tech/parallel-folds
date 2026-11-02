# -*- coding: utf-8 -*-
import os
import math
import json

from qgis.PyQt.QtCore import Qt, QSize, QPointF, QRectF, QMetaType
from qgis.PyQt.QtGui import QImage, QPainter, QColor, QPolygonF, QBrush, QPen, QFont
from qgis.PyQt.QtSvg import QSvgGenerator
from qgis.PyQt.QtWidgets import (
    QDialog, QWidget, QVBoxLayout, QHBoxLayout, QTableWidget, QTabWidget,
    QTableWidgetItem, QPushButton, QComboBox, QLabel, QFileDialog,
    QMessageBox, QCheckBox, QDoubleSpinBox, QGroupBox, QGridLayout,
    QFrame, QColorDialog
)

from qgis.core import (
    QgsProject, QgsRasterLayer, QgsVectorLayer, QgsFeature, QgsField,
    QgsGeometry, QgsPointXY, QgsWkbTypes, QgsRectangle,
    QgsCoordinateTransform, QgsCategorizedSymbolRenderer, QgsRendererCategory,
    QgsLineSymbol, QgsPalLayerSettings, QgsTextFormat, QgsVectorLayerSimpleLabeling,
    QgsUnitTypes, QgsSingleSymbolRenderer, QgsGraduatedSymbolRenderer,
    QgsVectorFileWriter, QgsMapLayer, QgsSymbol
)
from qgis.gui import QgsMapCanvas, QgsMapTool, QgsRubberBand, QgsVertexMarker

try:
    from qgis.core import Qgis
    JOIN_STYLE_ROUND = Qgis.JoinStyle.Round
except AttributeError:
    JOIN_STYLE_ROUND = QgsGeometry.JoinStyleRound

DEFAULT_COLORS = [
    QColor(220, 50, 50), QColor(50, 180, 50), QColor(50, 50, 220),
    QColor(230, 160, 10), QColor(160, 30, 230), QColor(30, 200, 230),
    QColor(230, 20, 150), QColor(120, 120, 120)
]


# =====================================================================
# STEREOPLOT CANVAS & DIALOG CLASS
# =====================================================================
class StereoplotCanvas(QWidget):
    """Advanced custom widget for stereoplot (Schmidt Equal-Area Projection)."""
    def __init__(self, bedding_data, parent=None):
        super().__init__(parent)
        self.bedding_data = bedding_data  # List: (feature_id, azimuth, dip)
        self.setMinimumSize(400, 400)
        
        self.show_poles = True
        self.show_planes = False
        self.show_contour = False
        self.show_beta = False
        
        self.beta_axis = None
        self._contour_image = None
        
        self.pole_vectors = []
        for fid, az, dip in self.bedding_data:
            dip_dir_rad = math.radians(az)
            dip_rad = math.radians(dip)
            x = -math.sin(dip_rad) * math.sin(dip_dir_rad)
            y = -math.sin(dip_rad) * math.cos(dip_dir_rad)
            z = -math.cos(dip_rad)
            self.pole_vectors.append((x, y, z))

    def update_contour_grid(self, r):
        if not self.show_contour or r <= 0:
            self._contour_image = None
            return

        side = int(r * 2)
        img = QImage(side, side, QImage.Format.Format_ARGB32)
        img.fill(Qt.GlobalColor.transparent)
        
        cos_max = 0.95
        max_count = 0
        counts = {}

        step = 4
        for dx in range(0, side, step):
            for dy in range(0, side, step):
                nx = (dx - r) / r
                ny = -(dy - r) / r
                r2 = nx*nx + ny*ny
                if r2 > 1.0:
                    continue
                
                sz = r2 - 1.0
                factor = math.sqrt(2.0 - r2)
                sx = nx * factor
                sy = ny * factor
                
                count = 0
                for px, py, pz in self.pole_vectors:
                    dot = sx*px + sy*py + sz*pz
                    if abs(dot) >= cos_max: 
                        count += 1
                
                if count > 0:
                    counts[(dx, dy)] = count
                    if count > max_count:
                        max_count = count

        if max_count > 0:
            painter = QPainter(img)
            for (dx, dy), count in counts.items():
                val = count / max_count
                color = QColor()
                color.setHsvF((1.0 - val) * 0.66, 0.9, 0.9, 0.5)
                painter.fillRect(dx, dy, step, step, QBrush(color))
            painter.end()

        self._contour_image = img

    def compute_beta_axis(self):
        if len(self.pole_vectors) < 2:
            self.beta_axis = None
            return

        m11, m12, m13 = 0.0, 0.0, 0.0
        m22, m23, m33 = 0.0, 0.0, 0.0

        for x, y, z in self.pole_vectors:
            m11 += x * x; m12 += x * y; m13 += x * z
            m22 += y * y; m23 += y * z; m33 += z * z

        A = [[m11, m12, m13], [m12, m22, m23], [m13, m23, m33]]
        V = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        
        for _ in range(50):
            row, col = 0, 1
            max_val = abs(A[0][1])
            if abs(A[0][2]) > max_val: row, col = 0, 2; max_val = abs(A[0][2])
            if abs(A[1][2]) > max_val: row, col = 1, 2; max_val = abs(A[1][2])

            if max_val < 1e-9: break

            phi = 0.5 * math.atan2(2.0 * A[row][col], A[row][row] - A[col][col])
            c, s = math.cos(phi), math.sin(phi)

            ar_r = c * c * A[row][row] + 2.0 * s * c * A[row][col] + s * s * A[col][col]
            ac_c = s * s * A[row][row] - 2.0 * s * c * A[row][col] + c * c * A[col][col]
            A[row][row], A[col][col] = ar_r, ac_c
            A[row][col] = A[col][row] = 0.0

            other = 3 - row - col
            a_ro = c * A[row][other] + s * A[col][other]
            a_co = -s * A[row][other] + c * A[col][other]
            A[row][other] = A[other][row] = a_ro
            A[col][other] = A[other][col] = a_co

            for k in range(3):
                v_kr = c * V[k][row] + s * V[k][col]
                v_kc = -s * V[k][row] + c * V[k][col]
                V[k][row], V[k][col] = v_kr, v_kc

        eigenvalues = [A[0][0], A[1][1], A[2][2]]
        min_idx = eigenvalues.index(min(eigenvalues))
        bx, by, bz = V[0][min_idx], V[1][min_idx], V[2][min_idx]

        if bz > 0: bx, by, bz = -bx, -by, -bz

        horiz_dist = math.hypot(bx, by)
        plunge = math.degrees(math.atan2(abs(bz), horiz_dist))
        trend = math.degrees(math.atan2(bx, by)) % 360.0

        self.beta_axis = (trend, plunge)

    def draw_stereonet(self, painter, width, height):
        cx, cy = width / 2, height / 2
        r = min(cx, cy) - 25

        if self.show_contour:
            if self._contour_image is None or self._contour_image.width() != int(r*2):
                self.update_contour_grid(r)
            if self._contour_image:
                painter.drawImage(QPointF(cx - r, cy - r), self._contour_image)

        pen_grid = QPen(QColor(160, 160, 160), 1, Qt.PenStyle.DashLine)
        pen_border = QPen(Qt.GlobalColor.black, 2)
        
        painter.setPen(pen_border)
        painter.drawEllipse(QPointF(cx, cy), r, r)
        painter.setPen(pen_grid)
        for grid_dip in [30, 60]:
            r_grid = r * math.sqrt(2) * math.sin(math.radians(grid_dip / 2.0))
            painter.drawEllipse(QPointF(cx, cy), r_grid, r_grid)
        painter.drawLine(QPointF(cx, cy - r), QPointF(cx, cy + r))
        painter.drawLine(QPointF(cx - r, cy), QPointF(cx + r, cy))

        def project_vector(vx, vy, vz):
            r_proj = r * math.sqrt(1.0 + vz)
            h_mag = math.hypot(vx, vy)
            if math.isclose(h_mag, 0.0): return cx, cy
            return cx + r_proj * (vx / h_mag), cy - r_proj * (vy / h_mag)

        if self.show_planes:
            painter.setPen(QPen(QColor(50, 50, 220, 140), 1.5))
            for fid, az, dip in self.bedding_data:
                dip_dir_rad, dip_rad = math.radians(az), math.radians(dip)
                poly = QPolygonF()
                for angle_deg in range(-90, 91, 2):
                    angle_rad = math.radians(angle_deg)
                    lx = math.sin(angle_rad)
                    ly = math.cos(angle_rad) * math.cos(dip_rad)
                    lz = -math.cos(angle_rad) * math.sin(dip_rad)
                    gx = lx * math.cos(dip_dir_rad) + ly * math.sin(dip_dir_rad)
                    gy = -lx * math.sin(dip_dir_rad) + ly * math.cos(dip_dir_rad)
                    px, py = project_vector(gx, gy, lz)
                    poly.append(QPointF(px, py))
                painter.drawPolyline(poly)

        if self.show_poles:
            painter.setPen(QPen(Qt.GlobalColor.red, 5))
            for x, y, z in self.pole_vectors:
                px, py = project_vector(x, y, z)
                painter.drawEllipse(QPointF(px, py), 3, 3)

        if self.show_beta and self.beta_axis:
            b_trend, b_plunge = self.beta_axis
            tr_rad, pl_rad = math.radians(b_trend), math.radians(b_plunge)
            bx = math.cos(pl_rad) * math.sin(tr_rad)
            by = math.cos(pl_rad) * math.cos(tr_rad)
            bz = -math.sin(pl_rad)
            px, py = project_vector(bx, by, bz)
            
            painter.setPen(QPen(Qt.GlobalColor.black, 2))
            painter.setBrush(QBrush(QColor(255, 235, 59)))
            painter.drawEllipse(QPointF(px, py), 7, 7)
            painter.setPen(QPen(Qt.GlobalColor.black, 4))
            painter.drawPoint(QPointF(px, py))

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.draw_stereonet(painter, self.width(), self.height())
        painter.end()


class StereoplotDialog(QDialog):
    def __init__(self, parent, bedding_data):
        super().__init__(parent)
        self.setWindowTitle("Structural Stereoplot (Schmidt Projection)")
        self.resize(500, 650)
        self.bedding_data = bedding_data  
        
        layout = QVBoxLayout(self)
        grid_buttons = QGridLayout()
        
        self.btn_poles = QPushButton("1. Plot Poles")
        self.btn_poles.setCheckable(True); self.btn_poles.setChecked(True)
        self.btn_poles.toggled.connect(self._toggle_poles)
        grid_buttons.addWidget(self.btn_poles, 0, 0)
        
        self.btn_planes = QPushButton("2. Plot Planes")
        self.btn_planes.setCheckable(True)
        self.btn_planes.toggled.connect(self._toggle_planes)
        grid_buttons.addWidget(self.btn_planes, 0, 1)
        
        self.btn_contour = QPushButton("3. Density Contour")
        self.btn_contour.setCheckable(True)
        self.btn_contour.toggled.connect(self._toggle_contour)
        grid_buttons.addWidget(self.btn_contour, 1, 0)
        
        self.btn_beta = QPushButton("4. Derive Beta Axis")
        self.btn_beta.setCheckable(True)
        self.btn_beta.toggled.connect(self._toggle_beta)
        grid_buttons.addWidget(self.btn_beta, 1, 1)
        
        layout.addLayout(grid_buttons)
        self.stereoplot_canvas = StereoplotCanvas(self.bedding_data, self)
        layout.addWidget(self.stereoplot_canvas, stretch=1)
        
        self.lbl_info = QLabel(f"Dataset: {len(self.bedding_data)} measured items.")
        self.lbl_info.setStyleSheet("font-weight: bold; font-size: 10pt; color: #004d40;")
        layout.addWidget(self.lbl_info)
        
        actions_layout = QHBoxLayout()
        self.btn_export_stereoplot_svg = QPushButton("Export Stereoplot SVG...")
        self.btn_export_stereoplot_svg.setStyleSheet("background-color: #37474f; color: white; font-weight: bold;")
        self.btn_export_stereoplot_svg.clicked.connect(self._export_stereoplot_to_svg)
        actions_layout.addWidget(self.btn_export_stereoplot_svg)
        
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        actions_layout.addWidget(btn_close)
        layout.addLayout(actions_layout)

    def _toggle_poles(self, checked):
        self.stereoplot_canvas.show_poles = checked; self.stereoplot_canvas.update()

    def _toggle_planes(self, checked):
        self.stereoplot_canvas.show_planes = checked; self.stereoplot_canvas.update()

    def _toggle_contour(self, checked):
        self.stereoplot_canvas.show_contour = checked; self.stereoplot_canvas._contour_image = None; self.stereoplot_canvas.update()

    def _toggle_beta(self, checked):
        self.stereoplot_canvas.show_beta = checked
        if checked:
            self.stereoplot_canvas.compute_beta_axis()
            res = self.stereoplot_canvas.beta_axis
            self.lbl_info.setText(f"Beta Axis: Trend {res[0]:.1f}° / Plunge {res[1]:.1f}°" if res else "Insufficient data.")
        else:
            self.lbl_info.setText(f"Dataset: {len(self.bedding_data)} measured items.")
        self.stereoplot_canvas.update()

    def _export_stereoplot_to_svg(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export as SVG", "stereoplot.svg", "Scalable Vector Graphics (*.svg)")
        if not path: return

        export_size = QSize(600, 600)
        generator = QSvgGenerator()
        generator.setFileName(path)
        generator.setSize(export_size)
        generator.setViewBox(QRectF(0, 0, export_size.width(), export_size.height()))

        painter = QPainter(generator)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(QRectF(0, 0, export_size.width(), export_size.height()), Qt.GlobalColor.white)
        self.stereoplot_canvas.draw_stereonet(painter, export_size.width(), export_size.height())
        painter.end()
        QMessageBox.information(self, "Exported", f"Stereoplot saved:\n{path}")


# =====================================================================
# MAIN PLUGIN CLASS: FOLDWINDOW
# =====================================================================
class FoldWindow(QDialog):
    
    def _sync_output_layer_to_store(self):
        """Synchronizes modified geometries directly respecting WebApp schema."""
        if not self._output_layer or not self.trace_layer:
            return

        sec_key = self._get_current_section_key()
        if not sec_key:
            return

        sec_layers_str, ok = QgsProject.instance().readEntry("ParallelFoldsPlugin", "section_layers", "")
        imported_data = json.loads(sec_layers_str) if (ok and sec_layers_str) else getattr(self, "_imported_section_layers", {})

        if "sectionStore" not in imported_data:
            imported_data["sectionStore"] = {}

        # Global Folds configuration
        folds_config = imported_data.get("foldsConfig", {
            "enabled": False,
            "layers": len(getattr(self, 'layers_thickness', [50.0, 50.0, 50.0])),
            "masterIdx": 0,
            "names": getattr(self, 'layer_names', ["Layer 1", "Layer 2", "Layer 3"]),
            "thicknesses": getattr(self, 'layers_thickness', [50.0, 50.0, 50.0]),
            "colors": [self.get_boundary_color(i).name() for i in range(len(getattr(self, 'layers_thickness', [50.0, 50.0, 50.0])))],
            "visible": [True] * len(getattr(self, 'layers_thickness', [50.0, 50.0, 50.0]))
        })

        updated_polylines = []
        for idx_feat, f in enumerate(self._output_layer.getFeatures()):
            geom = f.geometry()
            if geom.isEmpty(): 
                continue
            
            pts_json = [{"s": float(p.x()), "z": float(p.y())} for p in geom.vertices()]
            
            f_names = f.fields().names()
            l_idx = int(f["layer_idx"]) if "layer_idx" in f_names else (int(f["order"]) if "order" in f_names else idx_feat)
            col_hex = str(f["color"]) if "color" in f_names else "#1D9E75"
            
            master_idx = folds_config.get("masterIdx", 0)
            kind_str = "master" if l_idx == master_idx else "layer"

            updated_polylines.append({
                "id": idx_feat + 1,
                "kind": kind_str,
                "layerIdx": l_idx,
                "color": col_hex,
                "thickness": 2,
                "dashed": True,
                "pts": pts_json
            })

        # Save clean polylines per section
        imported_data["sectionStore"][sec_key] = {
            "polylines": updated_polylines
        }

        updated_json_str = json.dumps(imported_data)
        QgsProject.instance().writeEntry("ParallelFoldsPlugin", "section_layers", updated_json_str)
        self._imported_section_layers = imported_data

        if hasattr(self, '_compute_other_sections_intersections'):
            self._compute_other_sections_intersections()
        
        
    def _get_current_section_key(self):
        """Returns unique sec_key for the currently selected section line layer."""
        if not self.trace_layer:
            return None

        for feat in self.trace_layer.getFeatures():
            field_names = [f.name() for f in feat.fields()]
            if "sec_key" in field_names and feat.attribute("sec_key"):
                return str(feat.attribute("sec_key")).strip()
            
            geom = feat.geometry()
            if not geom.isEmpty():
                polyline = geom.asPolyline() if not geom.isMultipart() else geom.asMultiPolyline()[0]
                if len(polyline) >= 2:
                    pt_A, pt_B = polyline[0], polyline[-1]
                    return f"{pt_A.x():.2f}_{pt_A.y():.2f}_{pt_B.x():.2f}_{pt_B.y():.2f}"
            break
        return None
    
    
    def _compute_other_sections_intersections(self):
        """Calculates and displays intersection points of POLYLINES from other sections on current profile."""
        if not self.dem_layer or not self.trace_layer or not self._profile_layer:
            return

        section_features = list(self.trace_layer.getFeatures())
        if not section_features:
            return

        curr_geom = QgsGeometry(section_features[0].geometry())
        curr_geom.get().dropZValue(); curr_geom.get().dropMValue()

        xform_curr_to_dem = QgsCoordinateTransform(self.trace_layer.crs(), self.dem_layer.crs(), QgsProject.instance())
        curr_geom_dem = QgsGeometry(curr_geom)
        curr_geom_dem.transform(xform_curr_to_dem)

        sec_layers_str, ok = QgsProject.instance().readEntry("ParallelFoldsPlugin", "section_layers", "")
        imported_data = json.loads(sec_layers_str) if (ok and sec_layers_str) else getattr(self, "_imported_section_layers", {})
        sec_store = imported_data.get("sectionStore", {})

        intersections_found = []

        for layer_id, layer in QgsProject.instance().mapLayers().items():
            if not isinstance(layer, QgsVectorLayer) or layer.geometryType() != QgsWkbTypes.GeometryType.LineGeometry:
                continue
            
            if layer == self.trace_layer or layer in [self._profile_layer, self._output_layer, self._traces_layer, self._ghost_layer, self._strikes_layer, getattr(self, "_proj_intersection_lines_layer", None)]:
                continue

            xform_other_to_dem = QgsCoordinateTransform(layer.crs(), self.dem_layer.crs(), QgsProject.instance())

            for feat in layer.getFeatures():
                other_geom = QgsGeometry(feat.geometry())
                if other_geom.isEmpty(): continue
                other_geom.get().dropZValue(); other_geom.get().dropMValue()
                other_geom.transform(xform_other_to_dem)

                inter_geom = curr_geom_dem.intersection(other_geom)
                if inter_geom.isEmpty(): continue

                other_sec_key = None
                if "sec_key" in [f.name() for f in feat.fields()] and feat.attribute("sec_key"):
                    other_sec_key = str(feat.attribute("sec_key")).strip()
                else:
                    polyline_other = other_geom.asPolyline() if not other_geom.isMultipart() else other_geom.asMultiPolyline()[0]
                    if len(polyline_other) >= 2:
                        pA, pB = polyline_other[0], polyline_other[-1]
                        other_sec_key = f"{pA.x():.2f}_{pA.y():.2f}_{pB.x():.2f}_{pB.y():.2f}"

                other_polylines = sec_store.get(other_sec_key, {}).get("polylines", []) if other_sec_key else []

                for vertex in inter_geom.vertices():
                    pt_dem = QgsPointXY(vertex)
                    dist_s_curr = curr_geom_dem.lineLocatePoint(QgsGeometry.fromPointXY(pt_dem))
                    dist_s_other = other_geom.lineLocatePoint(QgsGeometry.fromPointXY(pt_dem))

                    if other_polylines:
                        for poly in other_polylines:
                            pts = poly.get("pts", poly.get("points", []))
                            if len(pts) < 2: continue
                            
                            z_interp = self._interpolate_z_at_s(pts, dist_s_other)
                            if z_interp is not None:
                                color_str = poly.get("color", "#00838f")
                                intersections_found.append((dist_s_curr, z_interp, color_str))
                    else:
                        dp_dem = self.dem_layer.dataProvider()
                        val, ok_sample = dp_dem.sample(pt_dem, 1)
                        qz = val if (ok_sample and not math.isnan(val)) else 0.0
                        intersections_found.append((dist_s_curr, qz, "#00838f"))

        if self._intersections_layer and self._intersections_layer.geometryType() != QgsWkbTypes.GeometryType.PointGeometry:
            self._remove_canvas_layer(self._intersections_layer)
            self._intersections_layer = None

        if not self._intersections_layer:
            crs_canvas = self.map_canvas.mapSettings().destinationCrs()
            self._intersections_layer = QgsVectorLayer(f"Point?crs={crs_canvas.authid()}", "Traces Intersection Points", "memory")
            provider = self._intersections_layer.dataProvider()
            provider.addAttributes([
                QgsField("id", QMetaType.Type.Int), QgsField("color", QMetaType.Type.QString)
            ])
            self._intersections_layer.updateFields()
            self._add_canvas_layer(self._intersections_layer)

        provider = self._intersections_layer.dataProvider()
        old_ids = [f.id() for f in self._intersections_layer.getFeatures()]
        if old_ids: 
            provider.deleteFeatures(old_ids)

        if intersections_found:
            feats_to_add = []
            categories = []
            used_colors = set()

            for dist_s, qz, col_str in intersections_found:
                feat = QgsFeature(self._intersections_layer.fields())
                feat.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(dist_s, qz)))
                feat.setAttributes([4, col_str])
                feats_to_add.append(feat)

                if col_str not in used_colors:
                    used_colors.add(col_str)
                    sym = QgsSymbol.defaultSymbol(self._intersections_layer.geometryType())
                    sym.setColor(QColor(col_str))
                    categories.append(QgsRendererCategory(col_str, sym, col_str))

            provider.addFeatures(feats_to_add)

            renderer = QgsCategorizedSymbolRenderer('color', categories)
            self._intersections_layer.setRenderer(renderer)

        self._intersections_layer.setLabelsEnabled(False)
        self._intersections_layer.updateExtents()
        self._intersections_layer.triggerRepaint()

        self._update_canvas_visibility()
        self.map_canvas.refresh()

    @staticmethod
    def _interpolate_z_at_s(pts, target_s):
        """Linearly interpolates Z elevation from a list of {s, z} points given target_s."""
        sorted_pts = sorted(pts, key=lambda p: p.get("s", 0.0))
        if target_s < sorted_pts[0]["s"] or target_s > sorted_pts[-1]["s"]:
            return None

        for i in range(len(sorted_pts) - 1):
            p1, p2 = sorted_pts[i], sorted_pts[i+1]
            s1, z1 = p1["s"], p1["z"]
            s2, z2 = p2["s"], p2["z"]

            if s1 <= target_s <= s2:
                if abs(s2 - s1) < 1e-6:
                    return z1
                t = (target_s - s1) / (s2 - s1)
                return z1 + t * (z2 - z1)

        return None
        
    def _on_tab_changed(self, index):
        if index != 1:
            if self.btn_start.isChecked(): self.btn_start.setChecked(False)
            if hasattr(self, 'btn_edit_saved') and self.btn_edit_saved.isChecked():
                self.btn_edit_saved.setChecked(False)
        
    def __init__(self, iface_obj):
        super().__init__(iface_obj.mainWindow())
        self.iface = iface_obj
        
        self.setWindowTitle("Parallel Folds & Structural Tool")
        self.resize(1150, 850)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)

        self.dem_layer = None
        self.trace_layer = None
        self.projected_structural_dips = []

        self.tool = None
        self.layer_names = ["Layer 1", "Layer 2", "Layer 3"]
        self.layers_thickness = [20.0, 20.0, 20.0]
        self.layers_colors = [DEFAULT_COLORS[i % len(DEFAULT_COLORS)] for i in range(16)]
        
        self.digitized_index = 0  
        self._canvas_layers = []
        self._output_layer = None
        self._ghost_layer = None  
        
        self._profile_layer = None       
        self._intersections_layer = None 
        self._strikes_layer = None       
        self._traces_layer = None        
        
        self._edit_tool = None
        self._select_proj_tool = None  

        self.profile_length = 1000.0
        self.profile_z_min = 0.0
        self.profile_z_max = 500.0
        self.topo_points = []  

        window_layout = QHBoxLayout()
        self.setLayout(window_layout)

        # 1. LEFT PANEL
        self.left_panel = QWidget()
        self.left_panel.setFixedWidth(340)
        left_layout = QVBoxLayout(self.left_panel)
        left_layout.setContentsMargins(0, 0, 0, 0)
        window_layout.addWidget(self.left_panel)

        self.main_tabs = QTabWidget()
        left_layout.addWidget(self.main_tabs, stretch=1)

        self.tab_profile_extraction = QWidget()
        self.tab_folding_digitization = QWidget()
        self.tab_import_export = QWidget()
        
        self.main_tabs.addTab(self.tab_profile_extraction, "1. Data Extraction & Profile")
        self.main_tabs.addTab(self.tab_folding_digitization, "2. Parallel Fold Digitization")
        self.main_tabs.addTab(self.tab_import_export, "3. Import / Export")

        self._setup_tab_profile()
        self._setup_tab_folding()
        self._setup_tab_import_export()

        # LAYER VISIBILITY CONTAINER
        self.visibility_container = QWidget()
        visibility_layout = QVBoxLayout(self.visibility_container)
        visibility_layout.setContentsMargins(15, 10, 15, 15)

        lbl_visibility = QLabel("Layers visibility")
        lbl_visibility.setStyleSheet("font-weight: bold; font-size: 11px;")
        visibility_layout.addWidget(lbl_visibility)

        self.chk_profile = QCheckBox("XY axes"); self.chk_profile.setChecked(True)
        self.chk_strikes = QCheckBox("Projected bedding"); self.chk_strikes.setChecked(True)
        self.chk_traces = QCheckBox("Projected traces"); self.chk_traces.setChecked(True)
        self.chk_intersections = QCheckBox("Traces intersection"); self.chk_intersections.setChecked(True)
        self.chk_parallel_folds = QCheckBox("Parallel Fold"); self.chk_parallel_folds.setChecked(True)

        visibility_layout.addWidget(self.chk_profile)
        visibility_layout.addWidget(self.chk_strikes)
        visibility_layout.addWidget(self.chk_traces)
        visibility_layout.addWidget(self.chk_intersections)
        visibility_layout.addWidget(self.chk_parallel_folds)

        self.chk_profile.toggled.connect(lambda c: self.toggle_layer_visibility("profile", c))
        self.chk_strikes.toggled.connect(lambda c: self.toggle_layer_visibility("strikes", c))
        self.chk_traces.toggled.connect(lambda c: self.toggle_layer_visibility("traces", c))
        self.chk_intersections.toggled.connect(lambda c: self.toggle_layer_visibility("intersections", c))
        self.chk_parallel_folds.toggled.connect(lambda c: self.toggle_layer_visibility("parallel_folds", c))

        left_layout.addWidget(self.visibility_container)

        # 2. QGIS CANVAS
        self.map_canvas = QgsMapCanvas()
        self.map_canvas.setCanvasColor(QColor(255, 255, 255))
        self.map_canvas.setDestinationCrs(QgsProject.instance().crs())
        self.map_canvas.setExtent(QgsRectangle(0, 0, 1000, 1000))
        self.map_canvas.enableAntiAliasing(True)
        window_layout.addWidget(self.map_canvas, stretch=1)

        QgsProject.instance().layersAdded.connect(self._populate_layers_combos)
        QgsProject.instance().layersRemoved.connect(self._populate_layers_combos)
        self.main_tabs.currentChanged.connect(self._on_tab_changed)

        self._refresh_table()
        self._refresh_combo()

    def toggle_layer_visibility(self, layer_key, visible):
        layer_map = {
            "profile": getattr(self, "_profile_layer", None),
            "strikes": getattr(self, "_strikes_layer", None),
            "traces": getattr(self, "_traces_layer", None),
            "intersections": getattr(self, "_intersections_layer", None),
            "parallel_folds": getattr(self, "_output_layer", None)
        }
        layer = layer_map.get(layer_key)
        if not layer: return

        root = QgsProject.instance().layerTreeRoot()
        layer_node = root.findLayer(layer.id())
        if layer_node:
            layer_node.setItemVisibilityChecked(visible)

        self._update_canvas_visibility()
        
    def _update_canvas_visibility(self):
        if not hasattr(self, "map_canvas") or not hasattr(self, "_canvas_layers"): return
        visible_layers = []
        for layer in self._canvas_layers:
            if not layer: continue
            if layer == getattr(self, "_profile_layer", None) and not self.chk_profile.isChecked(): continue
            if layer == getattr(self, "_strikes_layer", None) and not self.chk_strikes.isChecked(): continue
            if layer == getattr(self, "_traces_layer", None) and not self.chk_traces.isChecked(): continue
            if layer == getattr(self, "_intersections_layer", None) and not self.chk_intersections.isChecked(): continue
            if layer == getattr(self, "_proj_intersection_lines_layer", None) and not self.chk_intersections.isChecked(): continue
            if layer == getattr(self, "_output_layer", None) and not self.chk_parallel_folds.isChecked(): continue
            visible_layers.append(layer)

        self.map_canvas.setLayers(visible_layers)
        self.map_canvas.refresh()
    
    def _setup_tab_profile(self):
        layout = QVBoxLayout(self.tab_profile_extraction)
        self.sub_tabs = QTabWidget()
        layout.addWidget(self.sub_tabs)

        self.sub_tab_dem_trace = QWidget()
        self.sub_tab_2 = QWidget()
        self.sub_tab_project_traces = QWidget()  
        self.sub_tab_project_traces_four = QWidget()  
        
        self.sub_tabs.addTab(self.sub_tab_dem_trace, "DEM & Trace")
        self.sub_tabs.addTab(self.sub_tab_2, "Apparent Dip")
        self.sub_tabs.addTab(self.sub_tab_project_traces, "Trace Intersection")
        self.sub_tabs.addTab(self.sub_tab_project_traces_four, "Project Traces") 

        # -----------------------------------------------------------------
        # TAB 1: DEM & TRACE
        # -----------------------------------------------------------------
        sub_layout_1 = QVBoxLayout(self.sub_tab_dem_trace)
        sub_layout_1.addWidget(QLabel("Select Digital Elevation Model (DEM):"))
        self.combo_dem = QComboBox()
        sub_layout_1.addWidget(self.combo_dem)

        sub_layout_1.addWidget(QLabel("Select Section Line (Vector Layer):"))
        self.combo_trace = QComboBox()
        sub_layout_1.addWidget(self.combo_trace)

        sub_layout_1.addWidget(QLabel("Sampling Interval (Points count):"))
        self.spin_sampling = QDoubleSpinBox()
        self.spin_sampling.setRange(10.0, 5000.0); self.spin_sampling.setValue(200.0)
        sub_layout_1.addWidget(self.spin_sampling)

        self.btn_generate_profile = QPushButton("Extract & Generate Profile")
        self.btn_generate_profile.setStyleSheet("background-color: #2e7d32; color: white; font-weight: bold;")
        self.btn_generate_profile.clicked.connect(self._generate_topographic_profile)
        sub_layout_1.addWidget(self.btn_generate_profile)
        sub_layout_1.addStretch()

        # -----------------------------------------------------------------
        # TAB 2: APPARENT DIP (CORRETTO E COMPATIBILE QGIS 3.4)
        # -----------------------------------------------------------------
        layout_2 = QVBoxLayout(self.sub_tab_2)
        layout_2.setSpacing(6)
        layout_2.setContentsMargins(6, 6, 6, 6)

        layout_2.addWidget(QLabel("Select Structural Points Layer:"))
        self.combo_struct_layer = QComboBox()
        layout_2.addWidget(self.combo_struct_layer)
        self.combo_struct_layer.currentIndexChanged.connect(self._on_struct_layer_changed)

        layout_2.addWidget(QLabel("Layer Orientation Format:"))
        self.combo_orient_format = QComboBox()
        self.combo_orient_format.addItems(["Dip Direction / Dip", "Strike / Dip (Right-Hand Rule)"])
        layout_2.addWidget(self.combo_orient_format)

        # Layout affiancato per Dip Direction e Dip Field
        layout_dir_dip = QHBoxLayout()
        
        vbox_dir = QVBoxLayout()
        self.lbl_dir_field = QLabel("Dip Dir / Strike:")
        self.combo_dir_field = QComboBox()
        vbox_dir.addWidget(self.lbl_dir_field)
        vbox_dir.addWidget(self.combo_dir_field)
        layout_dir_dip.addLayout(vbox_dir)

        vbox_dip = QVBoxLayout()
        vbox_dip.addWidget(QLabel("Dip Field:"))
        self.combo_dip_field = QComboBox()
        vbox_dip.addWidget(self.combo_dip_field)
        layout_dir_dip.addLayout(vbox_dip)

        layout_2.addLayout(layout_dir_dip)

        # Gruppo Proiezione Personalizzata
        group_proj = QGroupBox("Custom Projection Vector")
        layout_proj = QGridLayout()
        layout_proj.addWidget(QLabel("Trend (0-360°):"), 0, 0)
        self.spin_proj_trend = QDoubleSpinBox()
        self.spin_proj_trend.setRange(0.0, 360.0); self.spin_proj_trend.setValue(90.0)
        layout_proj.addWidget(self.spin_proj_trend, 0, 1)

        layout_proj.addWidget(QLabel("Plunge (0-90°):"), 1, 0)
        self.spin_proj_plunge = QDoubleSpinBox()
        self.spin_proj_plunge.setRange(0.0, 90.0); self.spin_proj_plunge.setValue(0.0)
        layout_proj.addWidget(self.spin_proj_plunge, 1, 1)
        group_proj.setLayout(layout_proj)
        layout_2.addWidget(group_proj)

        # Parametri di ricerca e dimensione
        layout_h_params = QHBoxLayout()
        vbox_buf = QVBoxLayout()
        vbox_buf.addWidget(QLabel("Max Distance (m):"))
        self.spin_dip_buffer = QDoubleSpinBox()
        self.spin_dip_buffer.setRange(1.0, 100000.0); self.spin_dip_buffer.setValue(1000.0)
        vbox_buf.addWidget(self.spin_dip_buffer)
        layout_h_params.addLayout(vbox_buf)

        vbox_len = QVBoxLayout()
        vbox_len.addWidget(QLabel("Symbol Len (m):"))
        self.spin_dip_symbol_len = QDoubleSpinBox()
        self.spin_dip_symbol_len.setRange(5.0, 2000.0); self.spin_dip_symbol_len.setValue(50.0)
        vbox_len.addWidget(self.spin_dip_symbol_len)
        layout_h_params.addLayout(vbox_len)
        layout_2.addLayout(layout_h_params)

        # Pulsanti d'azione
        self.btn_project_dips = QPushButton("Project Apparent Dips")
        self.btn_project_dips.setStyleSheet("background-color: #ef6c00; color: white; font-weight: bold;")
        self.btn_project_dips.clicked.connect(self._compute_apparent_dips)
        layout_2.addWidget(self.btn_project_dips)
        
        self.btn_open_stereoplot = QPushButton("Open Stereoplot Window")
        self.btn_open_stereoplot.setStyleSheet("background-color: #00897b; color: white; font-weight: bold;")
        self.btn_open_stereoplot.clicked.connect(self._open_stereoplot)
        layout_2.addWidget(self.btn_open_stereoplot)
        
        layout_2.addStretch()

        # -----------------------------------------------------------------
        # TAB 3: TRACE INTERSECTION
        # -----------------------------------------------------------------
        layout_project_traces = QVBoxLayout(self.sub_tab_project_traces)
        layout_project_traces.addWidget(QLabel("Select Traces Layer (Lines/Polygons):"))
        self.combo_intersection_traces = QComboBox()
        layout_project_traces.addWidget(self.combo_intersection_traces)
        self.combo_intersection_traces.currentIndexChanged.connect(self._on_traces_layer_changed)

        layout_project_traces.addWidget(QLabel("Select Label Attribute Field:"))
        self.combo_trace_label_field = QComboBox()
        layout_project_traces.addWidget(self.combo_trace_label_field)

        self.btn_project_traces = QPushButton("Project Traces & Label Intersections")
        self.btn_project_traces.setStyleSheet("background-color: #1565c0; color: white; font-weight: bold;")
        self.btn_project_traces.clicked.connect(self._compute_intersections)
        layout_project_traces.addWidget(self.btn_project_traces)
        layout_project_traces.addStretch()

        # -----------------------------------------------------------------
        # TAB 4: PROJECT TRACES
        # -----------------------------------------------------------------
        layout_project_traces_four = QVBoxLayout(self.sub_tab_project_traces_four)
        layout_project_traces_four.addWidget(QLabel("Select Traces Layer (Lines/Polygons):"))
        self.combo_pt4_traces_layer = QComboBox()
        layout_project_traces_four.addWidget(self.combo_pt4_traces_layer)
        self.combo_pt4_traces_layer.currentIndexChanged.connect(self._on_pt4_traces_layer_changed)

        layout_project_traces_four.addWidget(QLabel("Select Label Attribute Field:"))
        self.combo_pt4_label_field = QComboBox()
        layout_project_traces_four.addWidget(self.combo_pt4_label_field)

        self.chk_pt4_selected_only = QCheckBox("Project selected features only in QGIS")
        layout_project_traces_four.addWidget(self.chk_pt4_selected_only)

        group_proj_four = QGroupBox("Projection Vector")
        layout_proj_four = QGridLayout()
        layout_proj_four.addWidget(QLabel("Trend (0-360°):"), 0, 0)
        self.spin_pt4_trend = QDoubleSpinBox()
        self.spin_pt4_trend.setRange(0.0, 360.0); self.spin_pt4_trend.setValue(90.0)
        layout_proj_four.addWidget(self.spin_pt4_trend, 0, 1)

        layout_proj_four.addWidget(QLabel("Plunge (0-90°):"), 1, 0)
        self.spin_pt4_plunge = QDoubleSpinBox()
        self.spin_pt4_plunge.setRange(0.0, 90.0); self.spin_pt4_plunge.setValue(0.0)
        layout_proj_four.addWidget(self.spin_pt4_plunge, 1, 1)
        group_proj_four.setLayout(layout_proj_four)
        layout_project_traces_four.addWidget(group_proj_four)

        layout_pt4_dist = QHBoxLayout()
        layout_pt4_dist.addWidget(QLabel("Max Distance (m):"))
        self.spin_pt4_buffer = QDoubleSpinBox()
        self.spin_pt4_buffer.setRange(1.0, 100000.0); self.spin_pt4_buffer.setValue(1000.0)
        layout_pt4_dist.addWidget(self.spin_pt4_buffer)
        layout_project_traces_four.addLayout(layout_pt4_dist)

        self.btn_pt4_project = QPushButton("Project Traces Nodes")
        self.btn_pt4_project.setStyleSheet("background-color: #6a1b9a; color: white; font-weight: bold;")
        self.btn_pt4_project.clicked.connect(self._compute_projected_traces_nodes)
        layout_project_traces_four.addWidget(self.btn_pt4_project)

        self.btn_pt4_select_single = QPushButton("Select & Delete Single Line")
        self.btn_pt4_select_single.setCheckable(True)
        self.btn_pt4_select_single.setStyleSheet("background-color: #0288d1; color: white; font-weight: bold;")
        self.btn_pt4_select_single.clicked.connect(self._toggle_select_proj_tool)
        layout_project_traces_four.addWidget(self.btn_pt4_select_single)

        self.btn_pt4_delete_selected = QPushButton("Delete Selected Line")
        self.btn_pt4_delete_selected.setStyleSheet("background-color: #e53935; color: white; font-weight: bold;")
        self.btn_pt4_delete_selected.clicked.connect(self._delete_selected_projected_line)
        layout_project_traces_four.addWidget(self.btn_pt4_delete_selected)

        self.btn_pt4_clear = QPushButton("Clear All Projected Traces")
        self.btn_pt4_clear.setStyleSheet("background-color: #c62828; color: white; font-weight: bold;")
        self.btn_pt4_clear.clicked.connect(self._clear_projected_traces_nodes)
        layout_project_traces_four.addWidget(self.btn_pt4_clear)

        layout_project_traces_four.addStretch()
        self._populate_layers_combos()

    def _setup_tab_folding(self):
        side = QVBoxLayout(self.tab_folding_digitization)

        self.chk_single_element = QCheckBox("Single layer only (Fault/Topography)")
        self.chk_single_element.toggled.connect(self._on_single_element_toggled)
        side.addWidget(self.chk_single_element)
        
        side.addWidget(QLabel("Layer thicknesses:"))
        self.table = QTableWidget(0, 2)
        self.table.setHorizontalHeaderLabels(["Layer", "Thickness"])
        self.table.setColumnWidth(0, 110); self.table.setColumnWidth(1, 120)
        side.addWidget(self.table)
        self.table.cellChanged.connect(self._on_thickness_changed)
        self.table.cellDoubleClicked.connect(self._change_layer_color)

        self.btn_row = QHBoxLayout()
        self.btn_add = QPushButton("+ layer"); self.btn_remove = QPushButton("- layer")
        self.btn_row.addWidget(self.btn_add); self.btn_row.addWidget(self.btn_remove)
        side.addLayout(self.btn_row)
        self.btn_add.clicked.connect(self._add_layer); self.btn_remove.clicked.connect(self._remove_layer)

        side.addWidget(QLabel("Boundary being digitized:"))
        self.combo_digitized = QComboBox()
        side.addWidget(self.combo_digitized)
        self.combo_digitized.currentIndexChanged.connect(self._on_digitized_changed)

        self.btn_start = QPushButton("Start digitalizing")
        self.btn_start.setCheckable(True)
        side.addWidget(self.btn_start)
        self.btn_start.toggled.connect(self._toggle_tool)

        self.btn_clear = QPushButton("Clear points")
        side.addWidget(self.btn_clear)
        self.btn_clear.clicked.connect(self._clear_points)

        self.btn_flip = QPushButton("Flip digitalizing direction (X)")
        side.addWidget(self.btn_flip)
        self.btn_flip.clicked.connect(self._flip_direction)

        self.btn_save = QPushButton("Save layers to project")
        side.addWidget(self.btn_save)
        self.btn_save.clicked.connect(self._save_layers)

        side.addWidget(QLabel("Saved layer editing:"))
        self.btn_edit_saved = QPushButton("Edit saved lines (Nodes)")
        self.btn_edit_saved.setCheckable(True)
        side.addWidget(self.btn_edit_saved)
        self.btn_edit_saved.toggled.connect(self._toggle_edit_tool)

        self.btn_load_selected = QPushButton("Load selected as master")
        side.addWidget(self.btn_load_selected)
        self.btn_load_selected.clicked.connect(self._load_selected_as_master)

        self.btn_delete_selected = QPushButton("Delete selected line(s)")
        side.addWidget(self.btn_delete_selected)
        self.btn_delete_selected.clicked.connect(self._delete_selected_lines)

        side.addStretch()

    def _setup_tab_import_export(self):
        """Setup for saving, loading, and exporting tools."""
        layout = QVBoxLayout(self.tab_import_export)

        # 1. GEOPACKAGE GROUP
        self.gpkg_group = QGroupBox("Legacy GeoPackage Storage")
        gpkg_layout = QVBoxLayout()
        
        self.btn_export_gpkg = QPushButton("Save All Layers to GeoPackage...")
        self.btn_export_gpkg.setVisible(False)
        self.btn_export_gpkg.setStyleSheet("background-color: #00796b; color: white; font-weight: bold;")
        self.btn_export_gpkg.clicked.connect(self._export_all_to_geopackage)
        gpkg_layout.addWidget(self.btn_export_gpkg)

        self.btn_import_gpkg = QPushButton("Load Layers from GeoPackage...")
        self.btn_import_gpkg.setStyleSheet("background-color: #004d40; color: white; font-weight: bold;")
        self.btn_import_gpkg.clicked.connect(self._import_from_geopackage)
        gpkg_layout.addWidget(self.btn_import_gpkg)
        
        self.gpkg_group.setLayout(gpkg_layout)
        layout.addWidget(self.gpkg_group)

        # 2. JSON GROUP (WEB-APP INTEGRATION)
        self.json_group = QGroupBox("GEOJSON Load/Save")
        json_layout = QVBoxLayout()

        self.btn_import_json = QPushButton("Load & Convert Field3Dnotes GeoJSON File")
        self.btn_import_json.setStyleSheet("background-color: #283593; color: white; font-weight: bold;")
        self.btn_import_json.clicked.connect(self._import_json_to_layers)
        json_layout.addWidget(self.btn_import_json)

        self.btn_export_json = QPushButton("Export Sections in GEOJSON")
        self.btn_export_json.setStyleSheet("background-color: #4a148c; color: white; font-weight: bold;")
        self.btn_export_json.clicked.connect(self._export_sections_to_json)
        json_layout.addWidget(self.btn_export_json)

        self.json_group.setLayout(json_layout)
        layout.addWidget(self.json_group)

        # SEPARATOR LINE
        linea_sep = QFrame()
        linea_sep.setFrameShape(QFrame.Shape.HLine)
        linea_sep.setFrameShadow(QFrame.Shadow.Sunken)
        layout.addWidget(linea_sep)

        # 3. GRAPHICAL EXPORT SECTION
        lbl_graphics = QLabel("Graphical Export")
        lbl_graphics.setStyleSheet("font-weight: bold;")
        layout.addWidget(lbl_graphics)

        self.btn_export_svg = QPushButton("Export to A4 SVG with Grid...")
        self.btn_export_svg.setStyleSheet("background-color: #37474f; color: white; font-weight: bold;")
        self.btn_export_svg.clicked.connect(self._export_svg)
        layout.addWidget(self.btn_export_svg)

        layout.addStretch()

    # =====================================================================
    # JSON / GEOJSON IMPORTER
    # =====================================================================
    def _import_json_to_layers(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self, 
            "Select JSON / GeoJSON File to Import", 
            "", 
            "JSON & GeoJSON Files (*.json *.geojson);;All Files (*)"
        )
        if not file_path: return

        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
        except Exception as e:
            QMessageBox.critical(self, "Loading Error", f"Unable to read JSON/GeoJSON file:\n{str(e)}")
            return

        crs_auth_id = "EPSG:3857"
        if isinstance(data, dict) and "crs" in data:
            crs_name = data["crs"].get("properties", {}).get("name", "")
            if "EPSG:" in crs_name.upper():
                crs_auth_id = "EPSG:" + crs_name.upper().split("EPSG:")[-1]

        features_list = data.get("features", [])

        structural_pts = []
        bestfit_centroids = []
        polylines = []
        draped_polygons = []
        sections = []

        for f in features_list:
            if not isinstance(f, dict): continue
            props = f.get("properties", {})
            l_type = props.get("layer_type", "")

            if l_type == "manual_orientation" or props.get("is_manual_spot", False):
                structural_pts.append(f)
            elif l_type == "orientation":
                polylines.append(f)
            elif l_type == "orientation_centroid":
                bestfit_centroids.append(f)
            elif l_type == "draped_polygon" or props.get("is_draped_polygon", False):
                draped_polygons.append(f)
            elif l_type == "simple_polyline" or props.get("is_simple_polyline", False):
                polylines.append(f)
            elif l_type == "geological_section" or props.get("is_geological_section", False):
                sections.append(f)

        if isinstance(data, dict) and "section_layers" in data:
            sec_layers_str = json.dumps(data["section_layers"])
            QgsProject.instance().writeEntry("ParallelFoldsPlugin", "section_layers", sec_layers_str)
            self._imported_section_layers = data["section_layers"]

        layers_created_count = 0

        # Layer 1: Structural Points
        if structural_pts:
            vl_struct = QgsVectorLayer(f"Point?crs={crs_auth_id}", "Structural Data Points", "memory")
            dp = vl_struct.dataProvider()
            dp.addAttributes([
                QgsField("id", QMetaType.Type.Int), QgsField("strike", QMetaType.Type.Double),
                QgsField("dip_dir", QMetaType.Type.Double), QgsField("dip", QMetaType.Type.Double)
            ])
            vl_struct.updateFields()
            q_feats = []
            for item in structural_pts:
                coords = item.get("geometry", {}).get("coordinates", [0, 0])
                feat = QgsFeature()
                feat.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(coords[0], coords[1])))
                p = item.get("properties", {})
                feat.setAttributes([
                    p.get("feature_id", 0), p.get("strike", 0.0),
                    p.get("dip_dir", 0.0), p.get("dip", 0.0)
                ])
                q_feats.append(feat)
            dp.addFeatures(q_feats); vl_struct.updateExtents()
            QgsProject.instance().addMapLayer(vl_struct); layers_created_count += 1

        # Layer 2: Centroids
        if bestfit_centroids:
            vl_cent = QgsVectorLayer(f"Point?crs={crs_auth_id}", "Best Fit Plane Centroids", "memory")
            dp = vl_cent.dataProvider()
            dp.addAttributes([
                QgsField("feature_id", QMetaType.Type.QString), QgsField("strike", QMetaType.Type.Double),
                QgsField("dip_dir", QMetaType.Type.Double), QgsField("dip", QMetaType.Type.Double)
            ])
            vl_cent.updateFields()
            q_feats = []
            for item in bestfit_centroids:
                coords = item.get("geometry", {}).get("coordinates", [0, 0])
                feat = QgsFeature()
                feat.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(coords[0], coords[1])))
                p = item.get("properties", {})
                feat.setAttributes([
                    str(p.get("feature_id", "")), p.get("strike", 0.0),
                    p.get("dip_dir", 0.0), p.get("dip", 0.0)
                ])
                q_feats.append(feat)
            dp.addFeatures(q_feats); vl_cent.updateExtents()
            QgsProject.instance().addMapLayer(vl_cent); layers_created_count += 1

        # Layer 3: Polylines & Section Traces
        if polylines:
            vl_lines = QgsVectorLayer(f"LineString?crs={crs_auth_id}", "Polylines", "memory")
            dp = vl_lines.dataProvider()
            dp.addAttributes([
                QgsField("id", QMetaType.Type.Int), QgsField("type", QMetaType.Type.QString)
            ])
            vl_lines.updateFields()
            q_feats = []
            for line_item in polylines:
                coords = line_item.get("geometry", {}).get("coordinates", [])
                if len(coords) < 2: continue
                q_pts = [QgsPointXY(c[0], c[1]) for c in coords]
                feat = QgsFeature()
                feat.setGeometry(QgsGeometry.fromPolylineXY(q_pts))
                p = line_item.get("properties", {})
                feat.setAttributes([p.get("feature_id", 0), p.get("layer_type", "")])
                q_feats.append(feat)
            dp.addFeatures(q_feats); vl_lines.updateExtents()
            QgsProject.instance().addMapLayer(vl_lines); layers_created_count += 1

        # Layer 4: Draped Polygons
        if draped_polygons:
            vl_poly = QgsVectorLayer(f"Polygon?crs={crs_auth_id}", "Draped Polygons", "memory")
            dp = vl_poly.dataProvider()
            dp.addAttributes([
                QgsField("id", QMetaType.Type.Int), QgsField("color", QMetaType.Type.QString)
            ])
            vl_poly.updateFields()
            q_feats = []
            for poly_item in draped_polygons:
                coords_rings = poly_item.get("geometry", {}).get("coordinates", [])
                if not coords_rings or len(coords_rings[0]) < 3: continue
                q_pts = [QgsPointXY(c[0], c[1]) for c in coords_rings[0]]
                feat = QgsFeature()
                feat.setGeometry(QgsGeometry.fromPolygonXY([q_pts]))
                p = poly_item.get("properties", {})
                feat.setAttributes([p.get("feature_id", 0), p.get("color", "#e60039")])
                q_feats.append(feat)
            dp.addFeatures(q_feats); vl_poly.updateExtents()
            
            first_color = draped_polygons[0].get("properties", {}).get("color", "#e60039")
            symbol = QgsSymbol.defaultSymbol(vl_poly.geometryType())
            symbol.setColor(QColor(first_color))
            symbol.setOpacity(0.5)
            vl_poly.setRenderer(QgsSingleSymbolRenderer(symbol))
            QgsProject.instance().addMapLayer(vl_poly); layers_created_count += 1

        # Layer 5: Geological Sections
        if sections:
            for idx, sec in enumerate(sections):
                p = sec.get("properties", {})
                sec_name = p.get("name") or f"Section_{p.get('feature_id', idx+1)}"
                vl_sec = QgsVectorLayer(f"LineString?crs={crs_auth_id}", f"Section - {sec_name}", "memory")
                dp_sec = vl_sec.dataProvider()
                dp_sec.addAttributes([
                    QgsField("feature_id", QMetaType.Type.Int),
                    QgsField("name", QMetaType.Type.QString),
                    QgsField("sec_key", QMetaType.Type.QString)
                ])
                vl_sec.updateFields()

                coords = sec.get("geometry", {}).get("coordinates", [])
                if len(coords) >= 2:
                    q_pts = [QgsPointXY(c[0], c[1]) for c in coords]
                    feat = QgsFeature()
                    feat.setGeometry(QgsGeometry.fromPolylineXY(q_pts))

                    A, B = p.get("A", {}), p.get("B", {})
                    sec_key = f"{A.get('E', 0):.2f}_{A.get('N', 0):.2f}_{B.get('E', 0):.2f}_{B.get('N', 0):.2f}"

                    feat.setAttributes([p.get("feature_id", idx + 1), sec_name, sec_key])
                    dp_sec.addFeatures([feat]); vl_sec.updateExtents()

                    line_symbol = QgsSymbol.defaultSymbol(vl_sec.geometryType())
                    line_symbol.setColor(QColor("#d32f2f")); line_symbol.setWidth(1.0)
                    vl_sec.setRenderer(QgsSingleSymbolRenderer(line_symbol))
                    QgsProject.instance().addMapLayer(vl_sec); layers_created_count += 1

        self._populate_layers_combos()
        QMessageBox.information(
            self, "Import Completed",
            f"Import completed successfully!\n\n"
            f"• Structural Data (Points): {len(structural_pts)}\n"
            f"• Centroids: {len(bestfit_centroids)}\n"
            f"• Polylines: {len(polylines)}\n"
            f"• Draped Polygons: {len(draped_polygons)}\n"
            f"• Geological Sections: {len(sections)}\n\n"
            f"Fold layers will be loaded automatically when selecting section lines."
        )

    # =====================================================================
    # EXPORT SECTIONS TO JSON (GLOBAL FOLDS CONFIG ONLY)
    # =====================================================================
    def _export_sections_to_json(self):
        """Exports GeoJSON with global foldsConfig and excluding empty sections."""
        
        sec_layers_str, ok = QgsProject.instance().readEntry("ParallelFoldsPlugin", "section_layers", "")
        current_data = json.loads(sec_layers_str) if (ok and sec_layers_str) else getattr(self, "_imported_section_layers", {})

        raw_thicknesses = getattr(self, 'layers_thickness', [50.0, 50.0, 50.0])
        default_colors = ["#1D9E75", "#e0a800", "#00e1ff", "#e83e8c", "#fd7e14"]
        
        folds_config = current_data.get("foldsConfig", {
            "enabled": False,
            "layers": len(raw_thicknesses),
            "masterIdx": 0,
            "names": getattr(self, 'layer_names', [f"Layer {i+1}" for i in range(len(raw_thicknesses))]),
            "thicknesses": [float(t) for t in raw_thicknesses],
            "colors": [self.get_boundary_color(i).name() if hasattr(self, 'get_boundary_color') else default_colors[i % len(default_colors)] for i in range(len(raw_thicknesses))],
            "visible": [True] * len(raw_thicknesses)
        })

        section_store = current_data.get("sectionStore", {})

        project_crs = QgsProject.instance().crs()
        crs_auth_id = project_crs.authid()
        epsg_code = crs_auth_id.split(":")[-1] if ":" in crs_auth_id else "3857"

        features_to_export = []
        clean_section_store = {}
        used_section_names = set()
        sec_keys_processed = set()
        store_keys = list(section_store.keys())
        
        global_feature_counter = 1

        for layer_id, layer in QgsProject.instance().mapLayers().items():
            if not isinstance(layer, QgsVectorLayer) or layer.geometryType() != QgsWkbTypes.GeometryType.LineGeometry:
                continue
            
            if layer in [self._profile_layer, self._output_layer, self._traces_layer, self._ghost_layer, self._strikes_layer, getattr(self, "_proj_intersection_lines_layer", None), getattr(self, "_intersections_layer", None)]:
                continue

            xform_to_project = QgsCoordinateTransform(layer.crs(), project_crs, QgsProject.instance()) if layer.crs() != project_crs else None

            for idx, feat in enumerate(layer.getFeatures()):
                geom = QgsGeometry(feat.geometry())
                if geom.isEmpty(): continue

                if xform_to_project:
                    geom.transform(xform_to_project)
                
                polyline = geom.asPolyline() if not geom.isMultipart() else geom.asMultiPolyline()[0]
                if len(polyline) < 2: continue

                pt_A, pt_B = polyline[0], polyline[-1]

                sec_key = f"{pt_A.x():.2f}_{pt_A.y():.2f}_{pt_B.x():.2f}_{pt_B.y():.2f}"
                proj_key = sec_key
                if sec_key not in section_store:
                    # keys saved by the plugin: layer 'sec_key' attribute or native-CRS endpoints
                    cands = []
                    if "sec_key" in feat.fields().names() and feat.attribute("sec_key"):
                        cands.append(str(feat.attribute("sec_key")).strip())
                    g_nat = feat.geometry()
                    if not g_nat.isEmpty():
                        pl_nat = g_nat.asPolyline() if not g_nat.isMultipart() else g_nat.asMultiPolyline()[0]
                        if len(pl_nat) >= 2:
                            cands.append(f"{pl_nat[0].x():.2f}_{pl_nat[0].y():.2f}_{pl_nat[-1].x():.2f}_{pl_nat[-1].y():.2f}")
                    sec_key = next((c for c in cands if c in section_store), sec_key)

                if sec_key not in section_store:
                    for s_k in store_keys:
                        parts = s_k.split("_")
                        if len(parts) == 4:
                            try:
                                ex, ey = float(parts[0]), float(parts[1])
                                if math.hypot(pt_A.x() - ex, pt_A.y() - ey) < 1.0:
                                    sec_key = s_k
                                    break
                            except ValueError: pass

                if sec_key in sec_keys_processed:
                    continue

                sec_data_entry = section_store.get(sec_key, {})
                polys = sec_data_entry.get("polylines", []) if isinstance(sec_data_entry, dict) else (sec_data_entry if isinstance(sec_data_entry, list) else [])

                if not polys:
                    continue

                sec_keys_processed.add(sec_key)

                clean_section_store[proj_key] = {
                    "polylines": polys
                }

                raw_name = feat.attribute("name") if "name" in [f.name() for f in feat.fields()] and feat.attribute("name") else layer.name()
                sec_name = str(raw_name).strip()
                if sec_name.lower().startswith("section - ") or sec_name.lower().startswith("sezione - "): 
                    sec_name = sec_name[10:].strip()

                final_name = sec_name
                counter = 1
                while final_name in used_section_names:
                    final_name = f"{sec_name}_{counter}"
                    counter += 1

                used_section_names.add(final_name)

                current_fid = global_feature_counter
                global_feature_counter += 1

                feature_dict = {
                    "type": "Feature",
                    "geometry": {
                        "type": "LineString",
                        "coordinates": [
                            [float(p.x()), float(p.y()), 0.0] for p in polyline
                        ]
                    },
                    "properties": {
                        "feature_id": current_fid,
                        "f_type": "Section",
                        "unit": "",
                        "set": "",
                        "color": "",
                        "object_type": "geological_section",
                        "layer_type": "geological_section",
                        "is_geological_section": True,
                        "name": final_name,
                        "A": { "E": float(pt_A.x()), "N": float(pt_A.y()), "Z": 0.0 },
                        "B": { "E": float(pt_B.x()), "N": float(pt_B.y()), "Z": 0.0 },
                        "section_data": None
                    }
                }
                features_to_export.append(feature_dict)

        if not features_to_export:
            QMessageBox.warning(self, "No Sections", "No active drawn sections available to export.")
            return

        path, _ = QFileDialog.getSaveFileName(
            self, 
            "Export Sections GeoJSON for Web-App", 
            "sections_export.geojson", 
            "GeoJSON Files (*.geojson *.json)"
        )
        if not path: return

        if not path.lower().endswith(".geojson") and not path.lower().endswith(".json"):
            path += ".geojson"

        export_geojson = {
            "type": "FeatureCollection",
            "crs": {
                "type": "name",
                "properties": { "name": f"urn:ogc:def:crs:EPSG::{epsg_code}" }
            },
            "section_layers": {
                "foldsConfig": folds_config,        
                "sectionStore": clean_section_store 
            },
            "features": features_to_export
        }

        try:
            with open(path, 'w', encoding='utf-8') as f:
                json.dump(export_geojson, f, indent=2, ensure_ascii=False)
            
            QMessageBox.information(
                self, 
                "Export Completed", 
                f"GeoJSON file saved successfully:\n{path}\n\n"
                f"• Valid exported sections: {len(features_to_export)}"
            )
        except Exception as e:
            QMessageBox.critical(self, "Save Error", f"Unable to save GeoJSON file:\n{str(e)}")
            
            
    # =====================================================================
    # FOLDS GENERATION AND PROJECTION LOGIC
    # =====================================================================
    def _generate_folds_for_selected_section(self, selected_layer):
        """Loads and extracts Parallel Folds for the active section."""
        if not selected_layer: 
            return

        if self._output_layer:
            old_ids = [f.id() for f in self._output_layer.getFeatures()]
            if old_ids:
                self._output_layer.dataProvider().deleteFeatures(old_ids)
                self._output_layer.triggerRepaint()

        sec_key = self._get_current_section_key()
        if not sec_key:
            self.map_canvas.refresh()
            return

        sec_layers_str, ok = QgsProject.instance().readEntry("ParallelFoldsPlugin", "section_layers", "")
        imported_data = json.loads(sec_layers_str) if (ok and sec_layers_str) else getattr(self, "_imported_section_layers", {})

        sec_store = imported_data.get("sectionStore", {})
        folds_config = imported_data.get("foldsConfig", {})

        self.layer_names = folds_config.get("names", getattr(self, 'layer_names', ["Layer 1", "Layer 2", "Layer 3"]))
        thicknesses = folds_config.get("thicknesses", getattr(self, 'layers_thickness', [50.0, 50.0, 50.0]))
        imported_colors = folds_config.get("colors", [])

        if imported_colors:
            for idx, hex_col in enumerate(imported_colors):
                if idx < len(self.layers_colors):
                    self.layers_colors[idx] = QColor(hex_col)

        if thicknesses:
            self.layers_thickness = [float(t) for t in thicknesses]
            self._refresh_table()
            self._refresh_combo()

        sec_data = sec_store.get(sec_key, {})
        folds_data = sec_data.get("polylines", []) if isinstance(sec_data, dict) else []

        if folds_data:
            crs_canvas = self.map_canvas.mapSettings().destinationCrs().authid()
            
            if not self._output_layer:
                self._output_layer = QgsVectorLayer(f"LineString?crs={crs_canvas}", "Parallel & Structural Layers", "memory")
                dp_fold = self._output_layer.dataProvider()
                dp_fold.addAttributes([
                    QgsField("layer_name", QMetaType.Type.QString),
                    QgsField("layer_idx", QMetaType.Type.Int),
                    QgsField("thickness", QMetaType.Type.Double),
                    QgsField("color", QMetaType.Type.QString)
                ])
                self._output_layer.updateFields()
                self._add_canvas_layer(self._output_layer)

            dp_fold = self._output_layer.dataProvider()
            feats_folds = []
            
            for p_item in folds_data:
                pts_raw = p_item.get("pts", p_item.get("points", []))
                if len(pts_raw) < 2: continue
                
                line_pts = [QgsPointXY(p["s"], p["z"]) for p in pts_raw if "s" in p and "z" in p]
                if len(line_pts) < 2: continue

                l_idx = p_item.get("layerIdx", p_item.get("layer_idx", 0))
                color_hex = p_item.get("color") or self.get_boundary_color(l_idx).name()
                l_name = self.layer_names[l_idx] if l_idx < len(self.layer_names) else f"Layer {l_idx+1}"
                l_thick = thicknesses[l_idx] if l_idx < len(thicknesses) else 50.0

                f = QgsFeature(self._output_layer.fields())
                f.setGeometry(QgsGeometry.fromPolylineXY(line_pts))
                
                field_names = [field.name() for field in self._output_layer.fields()]
                for f_name in field_names:
                    if f_name in ("layer_name", "label"):
                        f.setAttribute(f_name, l_name)
                    elif f_name in ("layer_idx", "order"):
                        f.setAttribute(f_name, l_idx)
                    elif f_name == "thickness":
                        f.setAttribute(f_name, l_thick)
                    elif f_name == "color":
                        f.setAttribute(f_name, color_hex)

                feats_folds.append(f)

            dp_fold.addFeatures(feats_folds)
            self._output_layer.updateExtents()
            self._apply_layer_symbology()

        self._update_canvas_visibility()
        self.map_canvas.refresh()

    def _open_stereoplot(self):
        if not self.trace_layer:
            QMessageBox.warning(self, "Missing trace", "Extract the topographic profile first to define the section line.")
            return

        struct_layer_id = self.combo_struct_layer.currentData()
        struct_layer = QgsProject.instance().mapLayer(struct_layer_id) if struct_layer_id else None
        dir_field = self.combo_dir_field.currentData()
        dip_field = self.combo_dip_field.currentData()

        if not struct_layer or not dir_field or not dip_field:
            QMessageBox.warning(self, "Missing Selection", "Verify that the structural layer and fields are selected.")
            return

        section_features = list(self.trace_layer.getFeatures())
        if not section_features: return

        section_geom = QgsGeometry(section_features[0].geometry())
        section_geom.get().dropZValue(); section_geom.get().dropMValue()

        xform_section_to_dem = QgsCoordinateTransform(self.trace_layer.crs(), self.dem_layer.crs(), QgsProject.instance())
        section_geom_dem_crs = QgsGeometry(section_geom)
        section_geom_dem_crs.transform(xform_section_to_dem)

        polyline = section_geom_dem_crs.asMultiPolyline()[0] if section_geom_dem_crs.isMultipart() else section_geom_dem_crs.asPolyline()
        if len(polyline) < 2: return

        buffer_dist = self.spin_dip_buffer.value()
        format_idx = self.combo_orient_format.currentIndex()
        section_length = section_geom_dem_crs.length()

        proj_trend = self.spin_proj_trend.value()
        proj_plunge = self.spin_proj_plunge.value()
        alpha_proj, beta_proj = math.radians(proj_trend), math.radians(proj_plunge)

        px = math.sin(alpha_proj) * math.cos(beta_proj)
        py = math.cos(alpha_proj) * math.cos(beta_proj)

        p_start, p_end = polyline[0], polyline[-1]
        A = p_end.y() - p_start.y()
        B = p_start.x() - p_end.x()
        C = p_end.x() * p_start.y() - p_start.x() * p_end.y()

        denom = A * px + B * py
        xform_struct_to_dem = QgsCoordinateTransform(struct_layer.crs(), self.dem_layer.crs(), QgsProject.instance())
        
        section_bbox = section_geom_dem_crs.boundingBox()
        section_bbox.grow(buffer_dist + 10.0)

        bedding_within_buffer = []

        for feature in struct_layer.getFeatures():
            geom = feature.geometry()
            if geom.isEmpty(): continue

            pt_dem_geom = QgsGeometry(geom)
            pt_dem_geom.transform(xform_struct_to_dem)
            pt_dem = pt_dem_geom.asPoint()

            if not section_bbox.contains(pt_dem): continue
            if section_geom_dem_crs.distance(pt_dem_geom) > buffer_dist: continue
            if math.isclose(denom, 0.0): continue

            t = -(A * pt_dem.x() + B * pt_dem.y() + C) / denom
            x_inter = pt_dem.x() + t * px
            y_inter = pt_dem.y() + t * py

            proj_geom_pt = QgsGeometry.fromPointXY(QgsPointXY(x_inter, y_inter))
            dist_local = section_geom_dem_crs.lineLocatePoint(proj_geom_pt)
            if dist_local < 0.0 or dist_local > section_length: continue

            try:
                raw_dir = float(feature[dir_field])
                dip = float(feature[dip_field])
                azimuth = raw_dir % 360.0 if format_idx == 0 else (raw_dir + 90.0) % 360.0
                bedding_within_buffer.append((feature.id(), azimuth, dip))
            except (ValueError, TypeError): continue

        if not bedding_within_buffer:
            QMessageBox.information(self, "No Data", "No bedding found within limits.")
            return

        self.stereoplot_window = StereoplotDialog(self, bedding_within_buffer)
        self.stereoplot_window.show()

    def _populate_layers_combos(self):
        if getattr(self, '_block_combo_updates', False): return

        prev_dem = self.combo_dem.currentData()
        prev_trace = self.combo_trace.currentData()
        prev_inter_trace = self.combo_intersection_traces.currentData()
        prev_pt4_trace = self.combo_pt4_traces_layer.currentData() if hasattr(self, 'combo_pt4_traces_layer') else None
        prev_struct = self.combo_struct_layer.currentData() if hasattr(self, 'combo_struct_layer') else None

        self.combo_dem.blockSignals(True)
        self.combo_trace.blockSignals(True)
        self.combo_intersection_traces.blockSignals(True)
        if hasattr(self, 'combo_pt4_traces_layer'): self.combo_pt4_traces_layer.blockSignals(True)
        if hasattr(self, 'combo_struct_layer'): self.combo_struct_layer.blockSignals(True)

        self.combo_dem.clear()
        self.combo_trace.clear()
        self.combo_intersection_traces.clear()
        if hasattr(self, 'combo_pt4_traces_layer'): self.combo_pt4_traces_layer.clear()
        if hasattr(self, 'combo_struct_layer'): self.combo_struct_layer.clear()

        excluded_keywords = {
            "topographic profile with grid", "projected bedding",
            "projected traces", "traces intersection",
            "profile base elements (ghost)", "parallel & structural layers"
        }

        for layer_id, layer in QgsProject.instance().mapLayers().items():
            layer_name_lower = layer.name().lower()
            if any(kw in layer_name_lower for kw in excluded_keywords): continue

            if isinstance(layer, QgsRasterLayer):
                if "background" in layer_name_lower or "calibrated" in layer_name_lower: continue
                self.combo_dem.addItem(layer.name(), layer_id)
            elif isinstance(layer, QgsVectorLayer):
                g_type = layer.geometryType()
                if g_type == 0 and hasattr(self, 'combo_struct_layer'):
                    self.combo_struct_layer.addItem(layer.name(), layer_id)
                elif g_type == 1:
                    self.combo_trace.addItem(layer.name(), layer_id)
                    self.combo_intersection_traces.addItem(layer.name(), layer_id)
                    if hasattr(self, 'combo_pt4_traces_layer'): self.combo_pt4_traces_layer.addItem(layer.name(), layer_id)
                elif g_type == 2:
                    self.combo_intersection_traces.addItem(layer.name(), layer_id)
                    if hasattr(self, 'combo_pt4_traces_layer'): self.combo_pt4_traces_layer.addItem(layer.name(), layer_id)

        def restore_selection(combo, prev_data):
            if prev_data:
                idx = combo.findData(prev_data)
                if idx >= 0: combo.setCurrentIndex(idx)

        restore_selection(self.combo_dem, prev_dem)
        restore_selection(self.combo_trace, prev_trace)
        restore_selection(self.combo_intersection_traces, prev_inter_trace)
        if hasattr(self, 'combo_pt4_traces_layer'): restore_selection(self.combo_pt4_traces_layer, prev_pt4_trace)
        if hasattr(self, 'combo_struct_layer'): restore_selection(self.combo_struct_layer, prev_struct)

        self.combo_dem.blockSignals(False)
        self.combo_trace.blockSignals(False)
        self.combo_intersection_traces.blockSignals(False)
        if hasattr(self, 'combo_pt4_traces_layer'): self.combo_pt4_traces_layer.blockSignals(False)
        if hasattr(self, 'combo_struct_layer'): self.combo_struct_layer.blockSignals(False)

        self._on_traces_layer_changed()
        self._on_pt4_traces_layer_changed()
        self._on_struct_layer_changed()
        
    def _on_traces_layer_changed(self):
        self.combo_trace_label_field.clear(); self.combo_trace_label_field.addItem("None", None)
        l_id = self.combo_intersection_traces.currentData()
        layer = QgsProject.instance().mapLayer(l_id) if l_id else None
        if layer and isinstance(layer, QgsVectorLayer):
            for f in layer.fields(): self.combo_trace_label_field.addItem(f.name(), f.name())

    def _on_pt4_traces_layer_changed(self):
        if not hasattr(self, 'combo_pt4_label_field'): return
        self.combo_pt4_label_field.clear(); self.combo_pt4_label_field.addItem("None (No Label)", None)
        l_id = self.combo_pt4_traces_layer.currentData()
        layer = QgsProject.instance().mapLayer(l_id) if l_id else None
        if layer and isinstance(layer, QgsVectorLayer):
            for f in layer.fields(): self.combo_pt4_label_field.addItem(f.name(), f.name())

    def _on_struct_layer_changed(self):
        if not hasattr(self, 'combo_struct_layer'): return
        prev_dir, prev_dip = self.combo_dir_field.currentData(), self.combo_dip_field.currentData()
        self.combo_dir_field.blockSignals(True); self.combo_dip_field.blockSignals(True)
        self.combo_dir_field.clear(); self.combo_dip_field.clear()

        l_id = self.combo_struct_layer.currentData()
        layer = QgsProject.instance().mapLayer(l_id) if l_id else None
        if layer and isinstance(layer, QgsVectorLayer):
            for f in layer.fields():
                if f.isNumeric():
                    self.combo_dir_field.addItem(f.name(), f.name())
                    self.combo_dip_field.addItem(f.name(), f.name())

        i = self.combo_dir_field.findData(prev_dir)
        if i >= 0: self.combo_dir_field.setCurrentIndex(i)
        j = self.combo_dip_field.findData(prev_dip)
        if j >= 0: self.combo_dip_field.setCurrentIndex(j)

        self.combo_dir_field.blockSignals(False); self.combo_dip_field.blockSignals(False)

    def _generate_topographic_profile(self):
        dem_layer_id, trace_layer_id = self.combo_dem.currentData(), self.combo_trace.currentData()
        if not dem_layer_id or not trace_layer_id:
            QMessageBox.warning(self, "Missing Layers", "Please select DEM and Trace layers.")
            return

        self.dem_layer = QgsProject.instance().mapLayer(dem_layer_id)
        self.trace_layer = QgsProject.instance().mapLayer(trace_layer_id)
        features = list(self.trace_layer.getFeatures())
       
        self.trace_layer = QgsProject.instance().mapLayer(trace_layer_id)
        if getattr(self, "_last_trace_id", None) != trace_layer_id:
            self._clear_projection_layers()
        self._last_trace_id = trace_layer_id
       
       
       
        if not features: return

        geom = QgsGeometry(features[0].geometry()); geom.get().dropZValue(); geom.get().dropMValue()
        xform = QgsCoordinateTransform(self.trace_layer.crs(), self.dem_layer.crs(), QgsProject.instance())
        geom_dem_crs = QgsGeometry(geom); geom_dem_crs.transform(xform)

        self.profile_length = geom_dem_crs.length()
        if self.profile_length <= 0: return

        num_samples = int(self.spin_sampling.value())
        self.topo_points = []
        dp = self.dem_layer.dataProvider()
        self.profile_z_min, self.profile_z_max = 999999.0, -999999.0

        for i in range(num_samples):
            dist_local = (i / (num_samples - 1)) * self.profile_length
            pt_geom = geom_dem_crs.interpolate(dist_local)
            if pt_geom.isEmpty(): continue
            pt_dem = pt_geom.asPoint()

            val, ok = dp.sample(pt_dem, 1)
            qz = val if (ok and not math.isnan(val)) else 0.0
            if qz < self.profile_z_min: self.profile_z_min = qz
            if qz > self.profile_z_max: self.profile_z_max = qz
            self.topo_points.append(QgsPointXY(dist_local, qz))

        if self._profile_layer:
            self._remove_canvas_layer(self._profile_layer); self._profile_layer = None

        crs_canvas = self.map_canvas.mapSettings().destinationCrs()
        self._profile_layer = QgsVectorLayer(f"LineString?crs={crs_canvas.authid()}", "Topographic Profile with Grid", "memory")
        provider = self._profile_layer.dataProvider()
        provider.addAttributes([
            QgsField("id", QMetaType.Type.Int), QgsField("label", QMetaType.Type.QString), QgsField("color", QMetaType.Type.QString)
        ])
        self._profile_layer.updateFields()

        feat_prof = QgsFeature(self._profile_layer.fields())
        feat_prof.setGeometry(QgsGeometry.fromPolylineXY(self.topo_points))
        feat_prof.setAttributes([1, "", ""])
        provider.addFeatures([feat_prof])

        self._update_grid_and_axes()
        self._update_profile_layer_renderer()

        label_settings = QgsPalLayerSettings()
        label_settings.fieldName = "label"
        text_format = QgsTextFormat(); text_format.setFont(QFont("Arial"))
        text_format.setSizeUnit(QgsUnitTypes.RenderUnit.RenderPixels); text_format.setSize(12.0); text_format.setColor(QColor("#333333"))
        label_settings.setFormat(text_format)
        label_settings.placement = QgsPalLayerSettings.Placement.Line
        label_settings.placementFlags = QgsPalLayerSettings.LinePlacementFlags.BelowLine
        label_settings.obstacle = False; label_settings.allowDegradedPlacement = True
        label_settings.priority = 10; label_settings.dist = 2.0 
        
        self._profile_layer.setLabeling(QgsVectorLayerSimpleLabeling(label_settings))
        self._profile_layer.setLabelsEnabled(True)

        self._add_canvas_layer(self._profile_layer, zoom_to_it=True)
        self._regenerate_ghost_layer()
        self._update_canvas_visibility()

        self._generate_folds_for_selected_section(self.trace_layer)
        self._compute_other_sections_intersections()
        
        
    def _update_grid_and_axes(self):
        if not self._profile_layer: return
        provider = self._profile_layer.dataProvider()
        
        old_ids = [f.id() for f in self._profile_layer.getFeatures() if f["id"] in (2, 3)]
        if old_ids: provider.deleteFeatures(old_ids)

        margin_y = (self.profile_z_max - self.profile_z_min) * 0.25 if (self.profile_z_max - self.profile_z_min) > 0 else 10.0
        axis_y_bottom = self.profile_z_min - margin_y
        axis_y_top = self.profile_z_max + margin_y

        len_x = self.profile_length
        len_y = axis_y_top - axis_y_bottom
        tick_size_comune = len_x / 50.0

        features_axes = []
        feat_ax = QgsFeature(self._profile_layer.fields())
        feat_ax.setGeometry(QgsGeometry.fromPolylineXY([QgsPointXY(0.0, axis_y_bottom), QgsPointXY(len_x, axis_y_bottom)]))
        feat_ax.setAttributes([2, "", ""]); features_axes.append(feat_ax)

        feat_ay = QgsFeature(self._profile_layer.fields())
        feat_ay.setGeometry(QgsGeometry.fromPolylineXY([QgsPointXY(0.0, axis_y_bottom), QgsPointXY(0.0, axis_y_top)]))
        feat_ay.setAttributes([2, "", ""]); features_axes.append(feat_ay)

        step_x = self._calculate_nice_step(len_x)
        curr_x = 0.0
        while curr_x <= len_x:
            feat_tick = QgsFeature(self._profile_layer.fields())
            feat_tick.setGeometry(QgsGeometry.fromPolylineXY([QgsPointXY(curr_x, axis_y_bottom), QgsPointXY(curr_x, axis_y_bottom - tick_size_comune)]))
            feat_tick.setAttributes([3, f"{curr_x:g}", ""]); features_axes.append(feat_tick)
            curr_x += step_x

        step_y = self._calculate_nice_step(len_y) / 0.5
        curr_y = math.ceil(axis_y_bottom / step_y) * step_y
        while curr_y <= axis_y_top:
            feat_tick_y = QgsFeature(self._profile_layer.fields())
            feat_tick_y.setGeometry(QgsGeometry.fromPolylineXY([QgsPointXY(0.0, curr_y), QgsPointXY(-tick_size_comune, curr_y)]))
            feat_tick_y.setAttributes([3, f"{curr_y:g}", ""]); features_axes.append(feat_tick_y)
            curr_y += step_y

        provider.addFeatures(features_axes)
        self._profile_layer.updateExtents(); self._profile_layer.updateFields(); self._profile_layer.triggerRepaint()
        
    def _compute_intersections(self):
        if not self.dem_layer or not self.trace_layer or not self._profile_layer:
            QMessageBox.warning(self, "No Profile", "Please extract profile first.")
            return

        traces_layer_id = self.combo_intersection_traces.currentData()
        traces_layer = QgsProject.instance().mapLayer(traces_layer_id) if traces_layer_id else None
        if not traces_layer: return

        section_features = list(self.trace_layer.getFeatures())
        if not section_features: return

        section_geom = QgsGeometry(section_features[0].geometry())
        section_geom.get().dropZValue(); section_geom.get().dropMValue()

        xform_section_to_dem = QgsCoordinateTransform(self.trace_layer.crs(), self.dem_layer.crs(), QgsProject.instance())
        section_geom_dem_crs = QgsGeometry(section_geom); section_geom_dem_crs.transform(xform_section_to_dem)

        xform_section_to_traces = QgsCoordinateTransform(self.trace_layer.crs(), traces_layer.crs(), QgsProject.instance())
        section_geom_traces_crs = QgsGeometry(section_geom); section_geom_traces_crs.transform(xform_section_to_traces)

        xform_traces_to_dem = QgsCoordinateTransform(traces_layer.crs(), self.dem_layer.crs(), QgsProject.instance())
        label_field = self.combo_trace_label_field.currentData()
        dp = self.dem_layer.dataProvider()

        intersections_found = []
        for feature in traces_layer.getFeatures():
            t_geom = feature.geometry()
            if t_geom.isEmpty(): continue
            t_geom.get().dropZValue(); t_geom.get().dropMValue()

            inter_geom = section_geom_traces_crs.intersection(t_geom)
            if inter_geom.isEmpty(): continue

            for vertex in inter_geom.vertices():
                pt_traces = QgsPointXY(vertex)
                pt_dem = xform_traces_to_dem.transform(pt_traces)
                
                dist_local = section_geom_dem_crs.lineLocatePoint(QgsGeometry.fromPointXY(pt_dem))
                val, ok = dp.sample(pt_dem, 1)
                qz = val if (ok and not math.isnan(val)) else 0.0

                label_val = str(feature[label_field]) if label_field and feature[label_field] is not None else f"ID {feature.id()}"
                intersections_found.append((dist_local, qz, label_val))

        proj_lines_layer = getattr(self, "_proj_intersection_lines_layer", None)
        if proj_lines_layer:
            self._remove_canvas_layer(proj_lines_layer)
            self._proj_intersection_lines_layer = None

        if intersections_found:
            crs_canvas = self.map_canvas.mapSettings().destinationCrs()
            self._proj_intersection_lines_layer = QgsVectorLayer(f"LineString?crs={crs_canvas.authid()}", "Projected Traces Intersections", "memory")
            provider = self._proj_intersection_lines_layer.dataProvider()
            provider.addAttributes([
                QgsField("id", QMetaType.Type.Int), QgsField("label", QMetaType.Type.QString), QgsField("color", QMetaType.Type.QString)
            ])
            self._proj_intersection_lines_layer.updateFields()

            margin_y = (self.profile_z_max - self.profile_z_min) * 0.05 if (self.profile_z_max - self.profile_z_min) > 0 else 10.0
            axis_y_bottom = self.profile_z_min - margin_y

            features_to_add = []
            for dist, qz, lbl in intersections_found:
                feat = QgsFeature(self._proj_intersection_lines_layer.fields())
                feat.setGeometry(QgsGeometry.fromPolylineXY([QgsPointXY(dist, axis_y_bottom), QgsPointXY(dist, qz)]))
                feat.setAttributes([4, lbl, ""]); features_to_add.append(feat)

            provider.addFeatures(features_to_add)
            
            sym = QgsLineSymbol.createSimple({'color': '#d32f2f', 'width': '1.0', 'line_style': 'dash'})
            self._proj_intersection_lines_layer.setRenderer(QgsSingleSymbolRenderer(sym))

            label_settings = QgsPalLayerSettings()
            label_settings.fieldName = "label"
            text_format = QgsTextFormat(); text_format.setFont(QFont("Arial")); text_format.setSizeUnit(QgsUnitTypes.RenderUnit.RenderPixels) 
            text_format.setSize(12.0); text_format.setColor(QColor("#d32f2f"))
            label_settings.setFormat(text_format)
            label_settings.placement = QgsPalLayerSettings.Placement.Line
            label_settings.placementFlags = QgsPalLayerSettings.LinePlacementFlags.AboveLine
            label_settings.obstacle = False; label_settings.allowDegradedPlacement = True; label_settings.priority = 8; label_settings.dist = 2.0 
            self._proj_intersection_lines_layer.setLabeling(QgsVectorLayerSimpleLabeling(label_settings))
            self._proj_intersection_lines_layer.setLabelsEnabled(True)

            self._proj_intersection_lines_layer.updateExtents()
            self._add_canvas_layer(self._proj_intersection_lines_layer)

        self._compute_other_sections_intersections()
        self._update_canvas_visibility()
        self.map_canvas.refresh()

    def _compute_apparent_dips(self):
        if not self.dem_layer or not self.trace_layer or not self._profile_layer: return

        struct_layer_id = self.combo_struct_layer.currentData()
        struct_layer = QgsProject.instance().mapLayer(struct_layer_id) if struct_layer_id else None
        dir_field, dip_field = self.combo_dir_field.currentData(), self.combo_dip_field.currentData()
        if not struct_layer or not dir_field or not dip_field: return

        section_features = list(self.trace_layer.getFeatures())
        if not section_features: return

        section_geom = QgsGeometry(section_features[0].geometry()); section_geom.get().dropZValue(); section_geom.get().dropMValue()
        xform_section_to_dem = QgsCoordinateTransform(self.trace_layer.crs(), self.dem_layer.crs(), QgsProject.instance())
        section_geom_dem_crs = QgsGeometry(section_geom); section_geom_dem_crs.transform(xform_section_to_dem)

        polyline = section_geom_dem_crs.asMultiPolyline()[0] if section_geom_dem_crs.isMultipart() else section_geom_dem_crs.asPolyline()
        if len(polyline) < 2: return

        dx_s = polyline[-1].x() - polyline[0].x()
        dy_s = polyline[-1].y() - polyline[0].y()
        section_azimuth = math.degrees(math.atan2(dx_s, dy_s)) % 360.0

        buffer_dist = self.spin_dip_buffer.value()
        symbol_len = self.spin_dip_symbol_len.value()
        format_idx = self.combo_orient_format.currentIndex()
        proj_trend, proj_plunge = self.spin_proj_trend.value(), self.spin_proj_plunge.value()

        xform_struct_to_dem = QgsCoordinateTransform(struct_layer.crs(), self.dem_layer.crs(), QgsProject.instance())
        dp = self.dem_layer.dataProvider()
        
        self.projected_structural_dips = []
        section_length = section_geom_dem_crs.length()
        section_bbox = section_geom_dem_crs.boundingBox(); section_bbox.grow(buffer_dist + 10.0)

        for feature in struct_layer.getFeatures():
            geom = feature.geometry()
            if geom.isEmpty(): continue

            pt_dem_geom = QgsGeometry(geom); pt_dem_geom.transform(xform_struct_to_dem)
            pt_dem = pt_dem_geom.asPoint()

            if not section_bbox.contains(pt_dem): continue
            if section_geom_dem_crs.distance(pt_dem_geom) > buffer_dist: continue

            try:
                raw_dir = float(feature[dir_field]); real_dip = float(feature[dip_field])
            except (ValueError, TypeError): continue

            dip_dir = raw_dir % 360.0 if format_idx == 0 else (raw_dir + 90.0) % 360.0
            angle_diff = math.radians(dip_dir - section_azimuth)
            apparent_dip_deg = math.degrees(math.atan(math.tan(math.radians(real_dip)) * math.cos(angle_diff)))

            alpha_proj, beta_proj = math.radians(proj_trend), math.radians(proj_plunge)
            px, py, pz = math.sin(alpha_proj) * math.cos(beta_proj), math.cos(alpha_proj) * math.cos(beta_proj), -math.sin(beta_proj)

            val_init, ok_init = dp.sample(pt_dem, 1)
            z_init = val_init if (ok_init and not math.isnan(val_init)) else 0.0

            p_start, p_end = polyline[0], polyline[-1]
            A = p_end.y() - p_start.y(); B = p_start.x() - p_end.x(); C = p_end.x() * p_start.y() - p_start.x() * p_end.y()

            denom = A * px + B * py
            if math.isclose(denom, 0.0): continue  

            t = -(A * pt_dem.x() + B * pt_dem.y() + C) / denom
            x_inter, y_inter, qz_projected = pt_dem.x() + t * px, pt_dem.y() + t * py, z_init + t * pz  

            proj_geom_pt = QgsGeometry.fromPointXY(QgsPointXY(x_inter, y_inter))
            if not section_geom_dem_crs.intersects(proj_geom_pt.buffer(0.001, 2)): continue 

            dist_local = section_geom_dem_crs.lineLocatePoint(proj_geom_pt)
            if dist_local < 0.0 or dist_local > section_length: continue

            self.projected_structural_dips.append((dist_local, qz_projected, apparent_dip_deg))

        if not self.projected_structural_dips: return

        if self._strikes_layer:
            self._remove_canvas_layer(self._strikes_layer); self._strikes_layer = None

        crs_canvas = self.map_canvas.mapSettings().destinationCrs()
        self._strikes_layer = QgsVectorLayer(f"LineString?crs={crs_canvas.authid()}", "Projected Bedding", "memory")
        provider = self._strikes_layer.dataProvider()
        provider.addAttributes([
            QgsField("id", QMetaType.Type.Int), QgsField("label", QMetaType.Type.QString), QgsField("color", QMetaType.Type.QString)
        ])
        self._strikes_layer.updateFields()

        features_to_add = []
        half_len = symbol_len / 2.0

        for dist, qz_projected, app_dip in self.projected_structural_dips:
            angle_rad = math.radians(app_dip)
            dx_p, dy_p = half_len * math.cos(angle_rad), -half_len * math.sin(angle_rad)

            feat_line = QgsFeature(self._strikes_layer.fields())
            feat_line.setGeometry(QgsGeometry.fromPolylineXY([QgsPointXY(dist - dx_p, qz_projected - dy_p), QgsPointXY(dist + dx_p, qz_projected + dy_p)]))
            feat_line.setAttributes([5, f"{abs(app_dip):.1f}°", ""]); features_to_add.append(feat_line)

            feat_dot = QgsFeature(self._strikes_layer.fields())
            feat_dot.setGeometry(QgsGeometry.fromPolylineXY([QgsPointXY(dist - symbol_len * 0.02, qz_projected), QgsPointXY(dist + symbol_len * 0.02, qz_projected)]))
            feat_dot.setAttributes([6, "", ""]); features_to_add.append(feat_dot)

        provider.addFeatures(features_to_add)

        categories = [
            QgsRendererCategory(5, QgsLineSymbol.createSimple({'color': '#2e7d32', 'width': '1.5'}), "Apparent Dip"),
            QgsRendererCategory(6, QgsLineSymbol.createSimple({'color': '#111111', 'width': '3.0'}), "Center Dot")
        ]
        self._strikes_layer.setRenderer(QgsCategorizedSymbolRenderer('id', categories))

        label_settings = QgsPalLayerSettings()
        label_settings.fieldName = "label"
        text_format = QgsTextFormat(); text_format.setFont(QFont("Arial")); text_format.setSizeUnit(QgsUnitTypes.RenderUnit.RenderPixels) 
        text_format.setSize(10.0); text_format.setColor(QColor("#2e7d32"))
        label_settings.setFormat(text_format)
        label_settings.placement = QgsPalLayerSettings.Placement.Line
        label_settings.placementFlags = QgsPalLayerSettings.LinePlacementFlags.AboveLine
        label_settings.obstacle = False; label_settings.dist = 4.0 
        self._strikes_layer.setLabeling(QgsVectorLayerSimpleLabeling(label_settings))
        self._strikes_layer.setLabelsEnabled(True)

        self._strikes_layer.updateExtents(); self._add_canvas_layer(self._strikes_layer); self._update_canvas_visibility()

    def _get_feature_color_safely(self, layer, feature):
        default_color = "#37474f"
        renderer = layer.renderer()
        if not renderer: return default_color
        try:
            if isinstance(renderer, QgsSingleSymbolRenderer):
                symbol = renderer.symbol()
                if symbol: return symbol.color().name()
            elif isinstance(renderer, QgsCategorizedSymbolRenderer):
                attr_name = renderer.classAttribute()
                val = feature[attr_name]
                for cat in renderer.categories():
                    if cat.value() == val:
                        sym = cat.symbol()
                        if sym: return sym.color().name()
            elif isinstance(renderer, QgsGraduatedSymbolRenderer):
                attr_name = renderer.classAttribute()
                try:
                    val = float(feature[attr_name])
                    for r in renderer.ranges():
                        if r.lowerValue() <= val <= r.upperValue():
                            sym = r.symbol()
                            if sym: return sym.color().name()
                except (ValueError, TypeError): pass
        except Exception: pass  # nosec B110
        return default_color

    def _compute_projected_traces_nodes(self):
        if not self.dem_layer or not self.trace_layer or not self._profile_layer: return

        pt4_layer_id = self.combo_pt4_traces_layer.currentData()
        pt4_layer = QgsProject.instance().mapLayer(pt4_layer_id) if pt4_layer_id else None
        if not pt4_layer: return

        features = list(pt4_layer.selectedFeatures()) if self.chk_pt4_selected_only.isChecked() else list(pt4_layer.getFeatures())
        section_features = list(self.trace_layer.getFeatures())
        if not section_features: return

        section_geom = QgsGeometry(section_features[0].geometry()); section_geom.get().dropZValue(); section_geom.get().dropMValue()
        xform_section_to_dem = QgsCoordinateTransform(self.trace_layer.crs(), self.dem_layer.crs(), QgsProject.instance())
        section_geom_dem_crs = QgsGeometry(section_geom); section_geom_dem_crs.transform(xform_section_to_dem)

        polyline = section_geom_dem_crs.asMultiPolyline()[0] if section_geom_dem_crs.isMultipart() else section_geom_dem_crs.asPolyline()
        if len(polyline) < 2: return

        buffer_dist, proj_trend, proj_plunge = self.spin_pt4_buffer.value(), self.spin_pt4_trend.value(), self.spin_pt4_plunge.value()
        xform_traces_to_dem = QgsCoordinateTransform(pt4_layer.crs(), self.dem_layer.crs(), QgsProject.instance())
        dp = self.dem_layer.dataProvider()
        section_length = section_geom_dem_crs.length()

        p_start, p_end = polyline[0], polyline[-1]
        A, B, C = p_end.y() - p_start.y(), p_start.x() - p_end.x(), p_end.x() * p_start.y() - p_start.x() * p_end.y()

        alpha_proj, beta_proj = math.radians(proj_trend), math.radians(proj_plunge)
        px, py, pz = math.sin(alpha_proj) * math.cos(beta_proj), math.cos(alpha_proj) * math.cos(beta_proj), -math.sin(beta_proj)

        denom = A * px + B * py
        if math.isclose(denom, 0.0): return

        chosen_label_field = self.combo_pt4_label_field.currentData()
        projected_lines_features = []

        crs_canvas = self.map_canvas.mapSettings().destinationCrs()
        temp_layer = QgsVectorLayer(f"LineString?crs={crs_canvas.authid()}", "Temp", "memory")
        temp_prov = temp_layer.dataProvider()
        temp_prov.addAttributes([
            QgsField("id", QMetaType.Type.Int), QgsField("label", QMetaType.Type.QString), QgsField("color", QMetaType.Type.QString)
        ])
        temp_layer.updateFields()

        for feature in features:
            feat_color_hex = self._get_feature_color_safely(pt4_layer, feature)
            lbl_val = str(feature[chosen_label_field]) if chosen_label_field and feature[chosen_label_field] is not None else ""

            t_geom = feature.geometry()
            if t_geom.isEmpty(): continue

            for part in t_geom.parts():
                vertices_list = list(part.vertices())
                if len(vertices_list) < 2: continue

                projected_points_in_part = []
                for i in range(len(vertices_list) - 1):
                    p1_raw, p2_raw = QgsPointXY(vertices_list[i]), QgsPointXY(vertices_list[i+1])

                    pt1_geom, pt2_geom = QgsGeometry.fromPointXY(p1_raw), QgsGeometry.fromPointXY(p2_raw)
                    pt1_geom.transform(xform_traces_to_dem); pt2_geom.transform(xform_traces_to_dem)
                    p1, p2 = pt1_geom.asPoint(), pt2_geom.asPoint()

                    if math.sqrt((p2.x() - p1.x())**2 + (p2.y() - p1.y())**2) < 1e-6: continue

                    for step in range(6):
                        t_ratio = step / 5.0
                        pt_curr = QgsPointXY(p1.x() + t_ratio * (p2.x() - p1.x()), p1.y() + t_ratio * (p2.y() - p1.y()))

                        val_init, ok_init = dp.sample(pt_curr, 1)
                        z_init = val_init if (ok_init and not math.isnan(val_init)) else 0.0

                        t_proj = -(A * pt_curr.x() + B * pt_curr.y() + C) / denom
                        if abs(t_proj) > buffer_dist: continue

                        proj_pt = QgsPointXY(pt_curr.x() + t_proj * px, pt_curr.y() + t_proj * py)
                        dist_local = section_geom_dem_crs.lineLocatePoint(QgsGeometry.fromPointXY(proj_pt))

                        if 0.0 <= dist_local <= section_length:
                            projected_points_in_part.append(QgsPointXY(dist_local, z_init + t_proj * pz))

                if len(projected_points_in_part) >= 2:
                    feat_line = QgsFeature(temp_layer.fields())
                    feat_line.setGeometry(QgsGeometry.fromPolylineXY(projected_points_in_part))
                    feat_line.setAttributes([7, lbl_val, feat_color_hex])
                    projected_lines_features.append(feat_line)

        if not projected_lines_features: return

        if self._traces_layer:
            self._remove_canvas_layer(self._traces_layer); self._traces_layer = None

        self._traces_layer = QgsVectorLayer(f"LineString?crs={crs_canvas.authid()}", "Projected Traces", "memory")
        provider = self._traces_layer.dataProvider()
        provider.addAttributes([
            QgsField("id", QMetaType.Type.Int), QgsField("label", QMetaType.Type.QString), QgsField("color", QMetaType.Type.QString)
        ])
        self._traces_layer.updateFields()
        provider.addFeatures(projected_lines_features)
        
        categories = []
        used_colors = {str(f["color"]) if f["color"] else "#37474f": True for f in self._traces_layer.getFeatures()}
        for c_hex in used_colors.keys():
            sym = QgsLineSymbol.createSimple({'color': c_hex, 'width': '0.8'})
            categories.append(QgsRendererCategory(c_hex, sym, f"Projected Trace ({c_hex})"))

        self._traces_layer.setRenderer(QgsCategorizedSymbolRenderer('color', categories))

        label_settings = QgsPalLayerSettings()
        label_settings.fieldName = "label"
        text_format = QgsTextFormat(); text_format.setFont(QFont("Arial")); text_format.setSizeUnit(QgsUnitTypes.RenderUnit.RenderPixels) 
        text_format.setSize(10.0); text_format.setColor(QColor("#333333"))
        label_settings.setFormat(text_format)
        label_settings.placement = QgsPalLayerSettings.Placement.Line
        label_settings.placementFlags = QgsPalLayerSettings.LinePlacementFlags.AboveLine
        label_settings.obstacle = False; label_settings.dist = 2.0 
        self._traces_layer.setLabeling(QgsVectorLayerSimpleLabeling(label_settings))
        self._traces_layer.setLabelsEnabled(True)

        self._traces_layer.updateExtents(); self._add_canvas_layer(self._traces_layer); self._update_canvas_visibility()

    def _update_profile_layer_renderer(self):
        if not self._profile_layer: return
        categories = [
            QgsRendererCategory("1", QgsLineSymbol.createSimple({'color': '#000000', 'width': '1.0'}), "Topography"),
            QgsRendererCategory("2", QgsLineSymbol.createSimple({'color': '#333333', 'width': '1.0'}), "Axes"),
            QgsRendererCategory("3", QgsLineSymbol.createSimple({'color': '#111111', 'width': '1.2'}), "Ticks")
        ]
        self._profile_layer.setRenderer(QgsCategorizedSymbolRenderer('to_string(id)', categories))
        self._profile_layer.triggerRepaint()

    def _toggle_select_proj_tool(self, checked):
        if checked:
            if self.btn_start.isChecked(): self.btn_start.setChecked(False)
            if self.btn_edit_saved.isChecked(): self.btn_edit_saved.setChecked(False)
            self._select_proj_tool = SelectProjectedLineTool(self.map_canvas, self)
            self.map_canvas.setMapTool(self._select_proj_tool)
            self.btn_pt4_select_single.setText("Stop Selecting Line")
        else:
            if self._select_proj_tool:
                self.map_canvas.unsetMapTool(self._select_proj_tool); self._select_proj_tool = None
            if self._traces_layer:
                self._traces_layer.removeSelection(); self.map_canvas.refresh()
            self.btn_pt4_select_single.setText("Select & Delete Single Projected Line")

    def _delete_selected_projected_line(self):
        if not self._traces_layer: return
        selected_ids = self._traces_layer.selectedFeatureIds()
        if not selected_ids: return

        to_delete = [fid for fid in selected_ids if self._traces_layer.getFeature(fid)["id"] == 7]
        if to_delete:
            self._traces_layer.dataProvider().deleteFeatures(to_delete)
            self._traces_layer.removeSelection(); self._traces_layer.updateExtents()
            self._traces_layer.triggerRepaint(); self.map_canvas.refresh()

    def _clear_projected_traces_nodes(self):
        if not self._traces_layer: return
        provider = self._traces_layer.dataProvider()
        old_ids = [f.id() for f in self._traces_layer.getFeatures() if f["id"] == 7]
        if old_ids:
            provider.deleteFeatures(old_ids); self._traces_layer.updateExtents()
            self._traces_layer.triggerRepaint(); self.map_canvas.refresh()

    def _regenerate_ghost_layer(self):
        if self._ghost_layer: self._remove_canvas_layer(self._ghost_layer); self._ghost_layer = None
        if not self.topo_points: return

        crs = self.map_canvas.mapSettings().destinationCrs().authid()
        self._ghost_layer = QgsVectorLayer(f"LineString?crs={crs}", "Profile Base Elements (Ghost)", "memory")
        provider = self._ghost_layer.dataProvider()
        provider.addAttributes([QgsField("element_type", QMetaType.Type.Int)])
        self._ghost_layer.updateFields()

        feat_topo = QgsFeature(self._ghost_layer.fields())
        feat_topo.setGeometry(QgsGeometry.fromPolylineXY(self.topo_points))
        feat_topo.setAttributes([1])
        provider.addFeatures([feat_topo]); self._ghost_layer.updateExtents()

        categories = [QgsRendererCategory(1, QgsLineSymbol.createSimple({'color': '#111111', 'width': '1.0'}), "Topography")]
        self._ghost_layer.setRenderer(QgsCategorizedSymbolRenderer('element_type', categories))
        self._add_canvas_layer(self._ghost_layer)

    def update_ghost_bounds_from_digitizing(self, min_y, max_y):
        c = False
        if min_y < self.profile_z_min: self.profile_z_min = min_y; c = True
        if max_y > self.profile_z_max: self.profile_z_max = max_y; c = True
        if c: self._regenerate_ghost_layer()

    def get_boundary_color(self, idx):
        return self.layers_colors[idx % len(self.layers_colors)]

    def _on_single_element_toggled(self, checked):
        self.table.setEnabled(not checked); self.btn_add.setEnabled(not checked)
        self.btn_remove.setEnabled(not checked); self.combo_digitized.setEnabled(not checked)
        self._live_update()

    def _add_canvas_layer(self, layer, zoom_to_it=False):
        QgsProject.instance().addMapLayer(layer, addToLegend=False)
        if layer not in self._canvas_layers: self._canvas_layers.insert(0, layer)
        self.map_canvas.setLayers(self._canvas_layers)
        if zoom_to_it: self.map_canvas.setExtent(layer.extent())
        self.map_canvas.refresh()

    def _remove_canvas_layer(self, layer):
        if layer in self._canvas_layers: self._canvas_layers.remove(layer)
        self.map_canvas.setLayers(self._canvas_layers)
        try: QgsProject.instance().removeMapLayer(layer.id())
        except Exception: pass  # nosec B110
        self.map_canvas.refresh()

    def _clear_projection_layers(self):
        """Removes layers derived from the previous section trace."""
        for attr in ("_strikes_layer", "_traces_layer",
                     "_intersections_layer", "_proj_intersection_lines_layer"):
            layer = getattr(self, attr, None)
            if layer:
                self._remove_canvas_layer(layer)
            setattr(self, attr, None)


    def _refresh_table(self):
        self.table.blockSignals(True)
        self.table.setRowCount(len(self.layers_thickness))
        
        layer_names = getattr(self, 'layer_names', [])
        
        for i, t in enumerate(self.layers_thickness):
            name_text = layer_names[i] if i < len(layer_names) else f"Layer {i+1}"
            
            item_name = QTableWidgetItem(name_text)
            item_name.setFlags(item_name.flags() | Qt.ItemFlag.ItemIsEditable)
            
            current_col = self.get_boundary_color(i)
            item_name.setBackground(QBrush(current_col))
            if current_col.lightness() < 128:
                item_name.setForeground(QBrush(QColor(255, 255, 255)))
            else:
                item_name.setForeground(QBrush(QColor(0, 0, 0)))

            self.table.setItem(i, 0, item_name)
            
            item_thick = QTableWidgetItem(str(t))
            item_thick.setFlags(item_thick.flags() | Qt.ItemFlag.ItemIsEditable)
            self.table.setItem(i, 1, item_thick)
            
        self.table.blockSignals(False)

    def _refresh_combo(self):
        self.combo_digitized.blockSignals(True); self.combo_digitized.clear()
        for i in range(len(self.layers_thickness) + 1):
            self.combo_digitized.addItem(f"Top Layer {i+1}" if i < len(self.layers_thickness) else f"Bottom Layer {i}")
        self.combo_digitized.setCurrentIndex(self.digitized_index)
        self.combo_digitized.blockSignals(False)

    def _add_layer(self):
        new_idx = len(self.layers_thickness) + 1
        self.layers_thickness.append(20.0)
        self.layer_names.append(f"Layer {new_idx}")
        self._refresh_table()
        self._refresh_combo()
        self._sync_folds_config()
        self._apply_layer_symbology()
        self._live_update()

    def _remove_layer(self):
        r = self.table.currentRow()
        if r < 0: r = len(self.layers_thickness) - 1
        if len(self.layers_thickness) > 1:
            del self.layers_thickness[r]
            if r < len(self.layer_names):
                del self.layer_names[r]
            self._refresh_table()
            self._refresh_combo()
            self._sync_folds_config()
            self._apply_layer_symbology()
            self._live_update()

    def _on_thickness_changed(self, row, col):
        if col == 0:
            new_name = self.table.item(row, 0).text()
            while len(self.layer_names) <= row:
                self.layer_names.append(f"Layer {len(self.layer_names)+1}")
            self.layer_names[row] = new_name
            self._sync_folds_config()
            self._apply_layer_symbology()
        elif col == 1:
            try:
                v = float(self.table.item(row, col).text().replace(",", "."))
                if v > 0: 
                    self.layers_thickness[row] = v
                    self._sync_folds_config()
            except (ValueError, AttributeError): pass
            self._live_update()

    def _change_layer_color(self, row, col):
        initial_color = self.get_boundary_color(row)
        color = QColorDialog.getColor(initial_color, self, f"Select color for {self.table.item(row, 0).text()}")
        if color.isValid():
            self.layers_colors[row % len(self.layers_colors)] = color
            self._refresh_table()
            self._sync_folds_config()
            self._apply_layer_symbology()
            self.map_canvas.refresh()

    def _sync_folds_config(self):
        sec_layers_str, ok = QgsProject.instance().readEntry("ParallelFoldsPlugin", "section_layers", "")
        imported_data = json.loads(sec_layers_str) if (ok and sec_layers_str) else getattr(self, "_imported_section_layers", {})

        if "foldsConfig" not in imported_data:
            imported_data["foldsConfig"] = {}

        imported_data["foldsConfig"]["names"] = self.layer_names[:len(self.layers_thickness)]
        imported_data["foldsConfig"]["thicknesses"] = self.layers_thickness
        imported_data["foldsConfig"]["colors"] = [self.get_boundary_color(i).name() for i in range(len(self.layers_thickness))]

        updated_str = json.dumps(imported_data)
        QgsProject.instance().writeEntry("ParallelFoldsPlugin", "section_layers", updated_str)
        self._imported_section_layers = imported_data

    def _live_update(self):
        if self.tool: self.tool.recalc_and_draw()

    def cumulative_offsets(self):
        if self.chk_single_element.isChecked(): return [0.0]
        nb = len(self.layers_thickness) + 1
        offsets = [0.0] * nb
        acc = 0.0
        for i in range(self.digitized_index + 1, nb):
            acc += self.layers_thickness[i - 1]; offsets[i] = -acc
        acc = 0.0
        for i in range(self.digitized_index - 1, -1, -1):
            acc += self.layers_thickness[i]; offsets[i] = acc
        return offsets

    def _on_digitized_changed(self, idx):
        if self.tool and len(self.tool.vertices) >= 2 and idx != self.digitized_index:
            m = next((g for (i, d, g) in self.tool.last_boundaries if i == idx), None)
            if m and not m.isEmpty():
                self.tool.set_vertices([QgsPointXY(p.x(), p.y()) for p in m.vertices()])
        self.digitized_index = idx
        self._live_update()

    def _flip_direction(self):
        if self.tool: self.tool.reverse_vertices()

    def _toggle_tool(self, checked):
        if checked:
            if self.btn_edit_saved.isChecked(): self.btn_edit_saved.setChecked(False)
            if hasattr(self, 'btn_pt4_select_single') and self.btn_pt4_select_single.isChecked():
                self.btn_pt4_select_single.setChecked(False); self._toggle_select_proj_tool(False)
            self.tool = ParallelFoldTool(self.map_canvas, self)
            self.map_canvas.setMapTool(self.tool)
            self.btn_start.setText("Stop digitalizing")
        else:
            if self.tool: self.tool.reset(); self.map_canvas.unsetMapTool(self.tool); self.tool = None
            self.btn_start.setText("Start digitalizing")

    def _apply_layer_symbology(self):
        if not self._output_layer: return
        
        field_names = [f.name() for f in self._output_layer.fields()]
        target_field = "layer_idx" if "layer_idx" in field_names else "order"
        
        layer_names = getattr(self, 'layer_names', [])
        
        categories = []
        for i in range(len(self.layers_thickness) + 1):
            lbl_name = layer_names[i] if i < len(layer_names) else f"Layer {i+1}"
            s = QgsLineSymbol.createSimple({'color': self.get_boundary_color(i).name(), 'width': '0.6'})
            categories.append(QgsRendererCategory(int(i), s, lbl_name))
            
        s_f = QgsLineSymbol.createSimple({'color': '#141414', 'width': '0.8'})
        categories.append(QgsRendererCategory(999, s_f, "Fault"))
        
        s_default = QgsLineSymbol.createSimple({'color': '#e53935', 'width': '0.6'})
        
        renderer = QgsCategorizedSymbolRenderer(target_field, categories)
        renderer.setSourceSymbol(s_default)
        
        self._output_layer.setRenderer(renderer)
        self._output_layer.triggerRepaint()

    def _save_layers(self):
        if not self.tool or len(self.tool.vertices) < 2: return
        self._block_combo_updates = True
        try:
            crs = self.map_canvas.mapSettings().destinationCrs().authid()
            
            if not self._output_layer:
                self._output_layer = QgsVectorLayer(f"LineString?crs={crs}", "Parallel & Structural Layers", "memory")
                dp = self._output_layer.dataProvider()
                dp.addAttributes([
                    QgsField("layer_name", QMetaType.Type.QString),
                    QgsField("layer_idx", QMetaType.Type.Int),
                    QgsField("thickness", QMetaType.Type.Double),
                    QgsField("color", QMetaType.Type.QString)
                ])
                self._output_layer.updateFields()
                self._add_canvas_layer(self._output_layer)

            field_names = [f.name() for f in self._output_layer.fields()]
            layer_names = getattr(self, 'layer_names', [])
            
            feats = []
            new_polylines_json = []

            for i, d, geom in self.tool.last_boundaries:
                f = QgsFeature(self._output_layer.fields())
                f.setGeometry(geom)
                
                idx_val = 999 if self.chk_single_element.isChecked() else int(i)
                color_val = str(self.get_boundary_color(i).name())
                
                if i < self.table.rowCount() and self.table.item(i, 0):
                    name_val = self.table.item(i, 0).text()
                elif i < len(layer_names):
                    name_val = layer_names[i]
                else:
                    name_val = f"Layer {i+1}"
                    
                thick_val = float(self.layers_thickness[i]) if i < len(self.layers_thickness) else 50.0

                for f_name in field_names:
                    if f_name in ("layer_name", "label"):
                        f.setAttribute(f_name, name_val)
                    elif f_name in ("layer_idx", "order"):
                        f.setAttribute(f_name, idx_val)
                    elif f_name == "thickness":
                        f.setAttribute(f_name, thick_val)
                    elif f_name == "color":
                        f.setAttribute(f_name, color_val)
                    elif f_name == "distance":
                        f.setAttribute(f_name, float(d) if d is not None else 0.0)

                feats.append(f)

                pts_json = [{"s": float(p.x()), "z": float(p.y())} for p in geom.vertices()]
                new_polylines_json.append({
                    "layerIdx": idx_val,
                    "color": color_val,
                    "pts": pts_json
                })

            self._output_layer.dataProvider().addFeatures(feats)
            self._output_layer.updateExtents()

            sec_key = self._get_current_section_key()
            if sec_key:
                sec_layers_str, ok = QgsProject.instance().readEntry("ParallelFoldsPlugin", "section_layers", "")
                imported_data = json.loads(sec_layers_str) if (ok and sec_layers_str) else getattr(self, "_imported_section_layers", {})

                if "sectionStore" not in imported_data:
                    imported_data["sectionStore"] = {}

                if sec_key not in imported_data["sectionStore"]:
                    imported_data["sectionStore"][sec_key] = {}

                imported_data["sectionStore"][sec_key].setdefault("polylines", []).extend(new_polylines_json)

                updated_json_str = json.dumps(imported_data)
                QgsProject.instance().writeEntry("ParallelFoldsPlugin", "section_layers", updated_json_str)
                self._imported_section_layers = imported_data

            self.tool.clear_points()
            self._apply_layer_symbology()
            self._update_canvas_visibility()
            self.map_canvas.refresh()
            self.btn_start.setChecked(False)
            
        finally:
            self._block_combo_updates = False
            self._populate_layers_combos()
            
    # GEOPACKAGE IMPORT / EXPORT
    def _export_all_to_geopackage(self):
        layers_to_export = {}
        if self._profile_layer and self._profile_layer.featureCount() > 0: layers_to_export["Profile_Layer"] = self._profile_layer
        if self._strikes_layer and self._strikes_layer.featureCount() > 0: layers_to_export["Strikes_Layer"] = self._strikes_layer
        if self._intersections_layer and self._intersections_layer.featureCount() > 0: layers_to_export["Intersections_Layer"] = self._intersections_layer
        if self._traces_layer and self._traces_layer.featureCount() > 0: layers_to_export["Traces_Layer"] = self._traces_layer
        if self._output_layer and self._output_layer.featureCount() > 0: layers_to_export["Parallel_Folds_Layer"] = self._output_layer

        if not layers_to_export: return
        path, _ = QFileDialog.getSaveFileName(self, "Save Workspace as GeoPackage", "section_workspace.gpkg", "GeoPackage (*.gpkg)")
        if not path: return

        if os.path.exists(path):
            try: os.remove(path)
            except Exception as e: pass  # nosec B110
        errors, is_first_layer = [], True
        transform_context = QgsProject.instance().transformContext()

        for name, layer in layers_to_export.items():
            options = QgsVectorFileWriter.SaveVectorOptions()
            options.fileEncoding = "UTF-8"; options.layerName = name; options.datasourceOptions = ["FORCE_FID=YES"]
            options.actionOnExistingFile = QgsVectorFileWriter.ActionOnExistingFile.CreateOrOverwriteFile if is_first_layer else QgsVectorFileWriter.ActionOnExistingFile.CreateOrOverwriteLayer
            is_first_layer = False

            err, err_msg, _, _ = QgsVectorFileWriter.writeAsVectorFormatV3(layer, path, transform_context, options)
            if err == QgsVectorFileWriter.WriterError.NoError:
                try:
                    from qgis.PyQt.QtXml import QDomDocument
                    doc = QDomDocument()
                    layer.exportNamedStyle(doc, QgsMapLayer.AllStyleCategories)
                    style_qml = doc.toString()
                    if style_qml:
                        layer.saveStyleToDatabase(name, "Default", "Workspace style", style_qml, True, "")
                except Exception: pass  # nosec B110
            else: errors.append(f"Layer '{name}' failed: {err_msg}")

        if not errors:
            QMessageBox.information(self, "Success", f"All workspace layers saved to:\n{path}")

    def _import_from_geopackage(self):
        path, _ = QFileDialog.getOpenFileName(self, "Load Workspace from GeoPackage", "", "GeoPackage (*.gpkg)")
        if not path: return

        gpkg_layers = {}
        import sqlite3
        try:
            conn = sqlite3.connect(path)
            cursor = conn.cursor()
            cursor.execute("SELECT table_name FROM gpkg_contents WHERE data_type='features'")
            for r in cursor.fetchall(): gpkg_layers[r[0]] = f"{path}|layername={r[0]}"
            conn.close()
        except Exception as e: return

        loaded_count = 0
        self._block_combo_updates = True

        try:
            if "Profile_Layer" in gpkg_layers:
                if self._profile_layer: self._remove_canvas_layer(self._profile_layer)
                self._profile_layer = QgsVectorLayer(gpkg_layers["Profile_Layer"], "Topographic Profile with Grid", "ogr")
                if self._profile_layer.isValid():
                    self.topo_points = []
                    for f in self._profile_layer.getFeatures():
                        if f["id"] == 1: self.topo_points = [QgsPointXY(pt) for pt in f.geometry().vertices()]; break
                    if self.topo_points:
                        self.profile_length = self.topo_points[-1].x() - self.topo_points[0].x()
                        self.profile_z_min = min(pt.y() for pt in self.topo_points)
                        self.profile_z_max = max(pt.y() for pt in self.topo_points)

                    self._update_profile_layer_renderer()
                    self._add_canvas_layer(self._profile_layer, zoom_to_it=True); self._regenerate_ghost_layer()
                    loaded_count += 1

            if "Strikes_Layer" in gpkg_layers:
                if self._strikes_layer: self._remove_canvas_layer(self._strikes_layer)
                self._strikes_layer = QgsVectorLayer(gpkg_layers["Strikes_Layer"], "Projected Bedding", "ogr")
                if self._strikes_layer.isValid():
                    self._add_canvas_layer(self._strikes_layer); loaded_count += 1

            if "Intersections_Layer" in gpkg_layers:
                if self._intersections_layer: self._remove_canvas_layer(self._intersections_layer)
                self._intersections_layer = QgsVectorLayer(gpkg_layers["Intersections_Layer"], "Traces Intersection Points", "ogr")
                if self._intersections_layer.isValid():
                    self._add_canvas_layer(self._intersections_layer); loaded_count += 1

            if "Traces_Layer" in gpkg_layers:
                if self._traces_layer: self._remove_canvas_layer(self._traces_layer)
                self._traces_layer = QgsVectorLayer(gpkg_layers["Traces_Layer"], "Projected Traces", "ogr")
                if self._traces_layer.isValid():
                    self._add_canvas_layer(self._traces_layer); loaded_count += 1

            if "Parallel_Folds_Layer" in gpkg_layers:
                if self._output_layer: self._remove_canvas_layer(self._output_layer)
                temp_ogr = QgsVectorLayer(gpkg_layers["Parallel_Folds_Layer"], "temp", "ogr")
                if temp_ogr.isValid():
                    crs = temp_ogr.crs().authid()
                    self._output_layer = QgsVectorLayer(f"LineString?crs={crs}", "Parallel & Structural Layers", "memory")
                    dp = self._output_layer.dataProvider()
                    dp.addAttributes([
                        QgsField("distance", QMetaType.Type.Double), QgsField("order", QMetaType.Type.Int), QgsField("color", QMetaType.Type.QString)
                    ])
                    self._output_layer.updateFields()
                    new_feats = []
                    for old_feat in temp_ogr.getFeatures():
                        nf = QgsFeature(self._output_layer.fields())
                        nf.setGeometry(old_feat.geometry())
                        nf["order"] = old_feat["order"] if "order" in old_feat.fields().names() else 0
                        nf["color"] = old_feat["color"] if "color" in old_feat.fields().names() else ""
                        new_feats.append(nf)
                    dp.addFeatures(new_feats); self._output_layer.updateExtents()
                    self._apply_layer_symbology(); self._add_canvas_layer(self._output_layer); loaded_count += 1

            self._update_canvas_visibility(); self.map_canvas.refresh()
        finally:
            self._block_combo_updates = False; self._populate_layers_combos()

    def _toggle_edit_tool(self, checked):
        if checked:
            if self.btn_start.isChecked(): self.btn_start.setChecked(False)
            if hasattr(self, 'btn_pt4_select_single') and self.btn_pt4_select_single.isChecked():
                self.btn_pt4_select_single.setChecked(False); self._toggle_select_proj_tool(False)
            self._edit_tool = AdvancedEditSavedLineTool(self.map_canvas, self)
            self.map_canvas.setMapTool(self._edit_tool)
        else:
            if self._edit_tool: 
                self._edit_tool.clear_markers(); self.map_canvas.unsetMapTool(self._edit_tool); self._edit_tool = None
            if self._output_layer: self._output_layer.removeSelection(); self._output_layer.triggerRepaint()
            self.map_canvas.refresh()
            
    def _load_selected_as_master(self):
        if not self._output_layer: return
        sel = list(self._output_layer.selectedFeatures())
        if len(sel) == 1:
            if not self.btn_start.isChecked(): self.btn_start.setChecked(True)
            self.tool.set_vertices([QgsPointXY(p.x(), p.y()) for p in sel[0].geometry().vertices()])

    def _delete_selected_lines(self):
        if not self._output_layer: return
        ids = self._output_layer.selectedFeatureIds()
        if not ids: return

        self._output_layer.dataProvider().deleteFeatures(ids)
        self._output_layer.updateExtents()
        self._output_layer.removeSelection()
        if self._edit_tool: self._edit_tool.clear_markers()

        if self.trace_layer:
            sec_key = None
            for feat in self.trace_layer.getFeatures():
                f_names = [field.name() for field in feat.fields()]
                if "sec_key" in f_names and feat.attribute("sec_key"):
                    sec_key = str(feat.attribute("sec_key")).strip()
                else:
                    geom_tr = feat.geometry()
                    if not geom_tr.isEmpty():
                        polyline = geom_tr.asPolyline() if not geom_tr.isMultipart() else geom_tr.asMultiPolyline()[0]
                        if len(polyline) >= 2:
                            pt_A, pt_B = polyline[0], polyline[-1]
                            sec_key = f"{pt_A.x():.2f}_{pt_A.y():.2f}_{pt_B.x():.2f}_{pt_B.y():.2f}"
                break

            if sec_key:
                sec_layers_str, ok = QgsProject.instance().readEntry("ParallelFoldsPlugin", "section_layers", "")
                imported_data = json.loads(sec_layers_str) if (ok and sec_layers_str) else getattr(self, "_imported_section_layers", {})
                
                sec_store = imported_data.get("sectionStore", {})
                if sec_key in sec_store:
                    rem_polylines = []
                    for f in self._output_layer.getFeatures():
                        if f.geometry().isEmpty(): continue
                        pts_json = [{"s": float(p.x()), "z": float(p.y())} for p in f.geometry().vertices()]
                        l_idx = f["layer_idx"] if "layer_idx" in f.fields().names() else 0
                        col_hex = f["color"] if "color" in f.fields().names() else "#e53935"
                        rem_polylines.append({"layerIdx": l_idx, "color": col_hex, "pts": pts_json})

                    sec_store[sec_key]["polylines"] = rem_polylines
                    imported_data["sectionStore"] = sec_store

                    updated_json_str = json.dumps(imported_data)
                    QgsProject.instance().writeEntry("ParallelFoldsPlugin", "section_layers", updated_json_str)
                    self._imported_section_layers = imported_data

        self.map_canvas.refresh()

    @staticmethod
    def _calculate_nice_step(delta):
        if delta is None or math.isnan(delta) or delta <= 0:
            return 10.0
        raw = delta / 6.0
        if raw <= 0 or math.isnan(raw):
            return 10.0
        log = math.floor(math.log10(raw))
        b = 10 ** log
        rel = raw / b
        nice = 1.0 if rel < 1.5 else 2.0 if rel < 3.5 else 5.0 if rel < 7.5 else 10.0
        return nice * b

    def _export_svg(self):
        data_to_export, texts_to_export, xs, ys = [], [], [], []

        if self._output_layer and self._output_layer.featureCount() > 0:
            for f in self._output_layer.getFeatures():
                if f.geometry().isEmpty(): continue
                pts = [QPointF(p.x(), p.y()) for p in f.geometry().vertices()]
                for p in f.geometry().vertices(): xs.append(p.x()); ys.append(p.y())
                data_to_export.append((pts, QColor(f["color"]), 0.5))

        if self._profile_layer and self._profile_layer.featureCount() > 0:
            for f in self._profile_layer.getFeatures():
                if f.geometry().isEmpty(): continue
                geom = f.geometry(); pts = [QPointF(p.x(), p.y()) for p in geom.vertices()]
                for p in geom.vertices(): xs.append(p.x()); ys.append(p.y())
                
                prof_id = f["id"]; color = QColor("#000000"); width = 0.4
                if prof_id == 1: width = 0.7
                elif prof_id == 2: color = QColor("#333333"); width = 0.3
                elif prof_id == 3:
                    color = QColor("#555555"); width = 0.2
                    fields_names = f.fields().names()
                    label = str(f["label"]) if "label" in fields_names and f["label"] is not None else ""
                    if label:
                        vertices = list(geom.vertices())
                        if len(vertices) >= 2:
                            p_text = vertices[1]
                            is_y_axis = (vertices[0].x() == vertices[1].x())
                            texts_to_export.append((label, p_text.x() + (-1.5 if is_y_axis else 0.0), p_text.y() + (0.0 if is_y_axis else -2.5), QColor("#333333"), 3.0, "right" if is_y_axis else "center"))
                data_to_export.append((pts, color, width))

        if self._intersections_layer and self._intersections_layer.featureCount() > 0:
            for f in self._intersections_layer.getFeatures():
                if f.geometry().isEmpty(): continue
                geom = f.geometry(); pts = [QPointF(p.x(), p.y()) for p in geom.vertices()]
                for p in geom.vertices(): xs.append(p.x()); ys.append(p.y())
                data_to_export.append((pts, QColor("#d32f2f"), 0.4))

                fields_names = f.fields().names()
                label = str(f["label"]) if "label" in fields_names and f["label"] is not None else ""
                if label:
                    vertices = list(geom.vertices())
                    if vertices:
                        top_pt = max(vertices, key=lambda p: p.y())
                        texts_to_export.append((label, top_pt.x(), top_pt.y() + 2.0, QColor("#d32f2f"), 3.5, "center"))

        if self._strikes_layer and self._strikes_layer.featureCount() > 0:
            for f in self._strikes_layer.getFeatures():
                if f.geometry().isEmpty(): continue
                geom = f.geometry(); pts = [QPointF(p.x(), p.y()) for p in geom.vertices()]
                for p in geom.vertices(): xs.append(p.x()); ys.append(p.y())
                col = QColor("#2e7d32") if f["id"] == 5 else QColor("#111111")
                width = 0.5 if f["id"] == 5 else 0.8
                data_to_export.append((pts, col, width))

                fields_names = f.fields().names()
                label = str(f["label"]) if "label" in fields_names and f["label"] is not None else ""
                
                
                if label and f["id"] == 5:
                    vertices = list(geom.vertices())
                    if vertices:
                        mid_x = sum(p.x() for p in vertices) / len(vertices)
                        mid_y = sum(p.y() for p in vertices) / len(vertices)
                        texts_to_export.append((label, mid_x, mid_y + 3.0, QColor("#2e7d32"), 2.8, "center"))

        if self._traces_layer and self._traces_layer.featureCount() > 0:
            for f in self._traces_layer.getFeatures():
                if f.geometry().isEmpty(): continue
                geom = f.geometry(); pts = [QPointF(p.x(), p.y()) for p in geom.vertices()]
                for p in geom.vertices(): xs.append(p.x()); ys.append(p.y())
                color = QColor(f["color"]) if f["color"] else QColor("#37474f")
                data_to_export.append((pts, color, 0.3))

                fields_names = f.fields().names()
                label = str(f["label"]) if "label" in fields_names and f["label"] is not None else ""

                if label:
                    vertices = list(geom.vertices())
                    if vertices:
                        mid_pt = vertices[len(vertices) // 2]
                        texts_to_export.append((label, mid_pt.x(), mid_pt.y() + 2.0, QColor("#333333"), 2.8, "center"))

        if not xs: return
        path, _ = QFileDialog.getSaveFileName(self, "Export SVG", "section_complete.svg", "SVG (*.svg)")
        if not path: return

        pw, ph, mar = 297.0, 210.0, 25.0
        x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
        sc = min((pw - 2*mar)/(x1-x0 or 1), (ph - 2*mar-10)/(y1-y0 or 1))
        ox = mar + ((pw-2*mar)-(x1-x0)*sc)/2
        oy = mar + ((ph-2*mar-10)-(y1-y0)*sc)/2 + 10

        def to_svg(x, y): return QPointF(ox + (x-x0)*sc, (ph-oy) - (y-y0)*sc)

        gen = QSvgGenerator()
        gen.setFileName(path); gen.setSize(QSize(int(pw*3.77), int(ph*3.77))); gen.setViewBox(QRectF(0, 0, pw, ph))
        
        p = QPainter()
        if p.begin(gen):
            p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
            p.fillRect(QRectF(0, 0, pw, ph), QBrush(Qt.GlobalColor.white))
            p.setPen(QPen(QColor(180, 180, 180), 0.15))
            p.drawRect(QRectF(to_svg(x0, y1), to_svg(x1, y0)))
            
            for pts, col, width in data_to_export:
                p.setPen(QPen(col, width)); p.drawPolyline(QPolygonF([to_svg(pt.x(), pt.y()) for pt in pts]))
            
            for text, x_map, y_map, col, size_mm, align in texts_to_export:
                pt_canvas = to_svg(x_map, y_map)
                font = QFont("Arial"); font.setPointSizeF(size_mm) 
                p.setFont(font); p.setPen(QPen(col, 0.1))
                fm = p.fontMetrics()
                tw, th = fm.horizontalAdvance(text), fm.height()
                text_x = pt_canvas.x() - (tw / 2.0) if align == "center" else pt_canvas.x() - tw - 1.0 if align == "right" else pt_canvas.x() + 1.0
                text_y = pt_canvas.y() + (th / 3.0) 
                p.drawText(QPointF(text_x, text_y), text)
            p.end()
            QMessageBox.information(self, "Success", "Profile section exported correctly to SVG.")
            
    def _clear_points(self):
        if self.tool: self.tool.clear_points()

    def closeEvent(self, event):
        if self.tool: self.map_canvas.unsetMapTool(self.tool)
        if self._edit_tool: self._edit_tool.clear_markers(); self.map_canvas.unsetMapTool(self._edit_tool)
        if self._select_proj_tool: self.map_canvas.unsetMapTool(self._select_proj_tool)
        for l in list(self._canvas_layers):
            try: QgsProject.instance().removeMapLayer(l.id())
            except: pass  # nosec B110
        super().closeEvent(event)


# =====================================================================
# INTERACTIVE MAP TOOLS
# =====================================================================
class ParallelFoldTool(QgsMapTool):
    def __init__(self, canvas, window):
        super().__init__(canvas)
        self.canvas, self.window = canvas, window
        self.vertices, self.vertex_markers, self.rubber_bands, self.last_boundaries = [], [], [], []
        self.dragged_index = None

    def clear_points(self):
        self.vertices = []
        for m in self.vertex_markers: self.canvas.scene().removeItem(m)
        self.vertex_markers = []
        for rb in self.rubber_bands: rb.reset(QgsWkbTypes.GeometryType.LineGeometry)
        self.rubber_bands = []

    def set_vertices(self, nv):
        self.clear_points(); self.vertices = list(nv)
        for pt in self.vertices:
            m = QgsVertexMarker(self.canvas); m.setCenter(pt)
            m.setIconType(QgsVertexMarker.ICON_CIRCLE); m.setIconSize(8); m.setPenWidth(2); m.setColor(QColor(55,138,221))
            self.vertex_markers.append(m)
        self.recalc_and_draw()

    def reverse_vertices(self):
        self.vertices.reverse(); self.vertex_markers.reverse(); self.recalc_and_draw()

    def reset(self): self.clear_points()

    def _hit_test(self, pos):
        for i, pt in enumerate(self.vertices):
            sc = self.toCanvasCoordinates(pt)
            if math.hypot(sc.x() - pos.x(), sc.y() - pos.y()) <= 12: return i
        return None

    def canvasPressEvent(self, e):
        if e.button() == Qt.MouseButton.RightButton: return
        idx = self._hit_test(e.pos())
        if idx is not None: self.dragged_index = idx; return
        
        pt = self.toMapCoordinates(e.pos()); pt_xy = QgsPointXY(pt)
        inserted = False
        
        if len(self.vertices) >= 2:
            geom_line = QgsGeometry.fromPolylineXY(self.vertices)
            tolerance = 12 * self.canvas.mapUnitsPerPixel()
            if geom_line.distance(QgsGeometry.fromPointXY(pt_xy)) <= tolerance:
                best_idx, min_dist = -1, float('inf')
                for i in range(len(self.vertices) - 1):
                    seg_geom = QgsGeometry.fromPolylineXY([self.vertices[i], self.vertices[i+1]])
                    dist = seg_geom.distance(QgsGeometry.fromPointXY(pt_xy))
                    if dist < min_dist: min_dist = dist; best_idx = i + 1
                
                if best_idx != -1:
                    self.vertices.insert(best_idx, pt_xy)
                    m = QgsVertexMarker(self.canvas); m.setCenter(pt_xy); m.setIconType(QgsVertexMarker.ICON_CIRCLE)
                    m.setIconSize(8); m.setColor(QColor(55,138,221)); self.vertex_markers.insert(best_idx, m)
                    inserted = True

        if not inserted:
            self.vertices.append(pt_xy)
            m = QgsVertexMarker(self.canvas); m.setCenter(pt_xy); m.setIconType(QgsVertexMarker.ICON_CIRCLE)
            m.setIconSize(8); m.setColor(QColor(55,138,221)); self.vertex_markers.append(m)
            
        self.recalc_and_draw()

    def canvasMoveEvent(self, e):
        if self.dragged_index is not None:
            pt = self.toMapCoordinates(e.pos())
            self.vertices[self.dragged_index] = QgsPointXY(pt)
            self.vertex_markers[self.dragged_index].setCenter(QgsPointXY(pt))
            self.recalc_and_draw()

    def canvasReleaseEvent(self, e): 
        self.dragged_index = None

    def canvasDoubleClickEvent(self, e):
        if self.feat and self.dragged_idx is not None and len(self.markers) > 2:
            m_to_remove = self.markers.pop(self.dragged_idx)
            self.canvas.scene().removeItem(m_to_remove)
            new_pts = [QgsPointXY(m.center()) for m in self.markers]
            self.window._output_layer.dataProvider().changeGeometryValues({self.feat.id(): QgsGeometry.fromPolylineXY(new_pts)})
            self.window._output_layer.triggerRepaint()
            self.canvas.refresh()
            self.dragged_marker, self.dragged_idx = None, None
            
            self.window._sync_output_layer_to_store()

    def recalc_and_draw(self):
        for rb in self.rubber_bands: rb.reset(QgsWkbTypes.GeometryType.LineGeometry)
        self.rubber_bands = []; self.last_boundaries = []
        if len(self.vertices) < 2: return
            
        self.window.update_ghost_bounds_from_digitizing(min(p.y() for p in self.vertices), max(p.y() for p in self.vertices))
        offsets = self.window.cumulative_offsets()
        lg = QgsGeometry.fromPolylineXY(self.vertices)
        sign = 1.0 if self.vertices[-1].x() >= self.vertices[0].x() else -1.0
        
        for i, d in enumerate(offsets):
            ad = d * sign
            geom = QgsGeometry(lg) if abs(ad) < 1e-9 else lg.offsetCurve(ad, 8, JOIN_STYLE_ROUND, 2.0)
            if not geom or geom.isEmpty(): continue
            self.last_boundaries.append((i, d, geom))
            rb = QgsRubberBand(self.canvas, QgsWkbTypes.GeometryType.LineGeometry)
            rb.setToGeometry(geom, None); rb.setColor(self.window.get_boundary_color(i))
            rb.setWidth(4 if i == self.window.digitized_index else 2)
            self.rubber_bands.append(rb)


class AdvancedEditSavedLineTool(QgsMapTool):
    def __init__(self, canvas, window):
        super().__init__(canvas)
        self.canvas, self.window = canvas, window
        self.feat, self.markers, self.dragged_marker, self.dragged_idx = None, [], None, None

    def clear_markers(self):
        for m in self.markers:
            if m: self.canvas.scene().removeItem(m)
        self.markers = []

    def deactivate(self):
        self.window._sync_output_layer_to_store()
        self.clear_markers()
        if self.window._output_layer:
            self.window._output_layer.removeSelection()
            self.window._output_layer.triggerRepaint()
        self.canvas.refresh()
        super().deactivate()

    def canvasPressEvent(self, e):
        if self.window.main_tabs.currentIndex() != 1: return
        layer = self.window._output_layer
        if not layer: return
        pt = self.toMapCoordinates(e.pos())
        
        for i, m in enumerate(self.markers):
            if math.hypot(m.center().x() - pt.x(), m.center().y() - pt.y()) < 15 * self.canvas.mapUnitsPerPixel():
                self.dragged_marker, self.dragged_idx = m, i
                return

        cg = QgsGeometry.fromPointXY(pt)
        bf, bd = None, 12 * self.canvas.mapUnitsPerPixel()
        for f in layer.getFeatures():
            d = f.geometry().distance(cg)
            if d <= bd: bd, bf = d, f
        
        if bf:
            layer.selectByIds([bf.id()]); self.feat = bf; self.clear_markers()
            for p in bf.geometry().vertices():
                m = QgsVertexMarker(self.canvas); m.setCenter(QgsPointXY(p))
                m.setIconType(QgsVertexMarker.ICON_BOX); m.setColor(QColor(230, 20, 150))
                self.markers.append(m)
        else:
            layer.removeSelection(); self.feat = None; self.clear_markers()
        
        layer.triggerRepaint(); self.canvas.refresh()

    def canvasMoveEvent(self, e):
        if self.dragged_marker:
            pt = self.toMapCoordinates(e.pos())
            self.dragged_marker.setCenter(QgsPointXY(pt))

    def canvasReleaseEvent(self, e):
        if self.dragged_marker and self.feat:
            new_pts = [QgsPointXY(m.center()) for m in self.markers]
            
            self.window._output_layer.dataProvider().changeGeometryValues({self.feat.id(): QgsGeometry.fromPolylineXY(new_pts)})
            self.window._output_layer.updateExtents()
            self.window._output_layer.triggerRepaint()
            self.canvas.refresh()
            self.dragged_marker = None
            
            self.window._sync_output_layer_to_store()

    def canvasDoubleClickEvent(self, e):
        if self.feat and self.dragged_idx is not None and len(self.markers) > 2:
            m_to_remove = self.markers.pop(self.dragged_idx)
            self.canvas.scene().removeItem(m_to_remove)
            new_pts = [QgsPointXY(m.center()) for m in self.markers]
            
            self.window._output_layer.dataProvider().changeGeometryValues({self.feat.id(): QgsGeometry.fromPolylineXY(new_pts)})
            self.window._output_layer.updateExtents()
            self.window._output_layer.triggerRepaint()
            self.canvas.refresh()
            self.dragged_marker, self.dragged_idx = None, None
            
            self.window._sync_output_layer_to_store()

class SelectProjectedLineTool(QgsMapTool):
    def __init__(self, canvas, window):
        super().__init__(canvas)
        self.canvas, self.window = canvas, window

    def canvasPressEvent(self, e):
        layer = self.window._traces_layer
        if not layer: return

        pt = self.toMapCoordinates(e.pos())
        cg = QgsGeometry.fromPointXY(pt)
        search_distance = 10 * self.canvas.mapUnitsPerPixel()
        best_feat, min_dist = None, search_distance

        for f in layer.getFeatures():
            if f["id"] != 7: continue
            dist = f.geometry().distance(cg)
            if dist <= min_dist: min_dist, best_feat = dist, f

        if best_feat: layer.selectByIds([best_feat.id()])
        else: layer.removeSelection()
            
        self.canvas.refresh()

    def keyPressEvent(self, e):
        if e.key() in (Qt.Key.Key_Delete, Qt.Key.Key_Backspace):
            self.window._delete_selected_projected_line()
