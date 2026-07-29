# -*- coding: utf-8 -*-
import os
import math

from qgis.PyQt.QtCore import Qt, QSize, QPointF, QRectF, QMetaType
from qgis.PyQt.QtGui import QImage, QPainter, QColor, QPolygonF, QBrush, QPen, QFont
from qgis.PyQt.QtSvg import QSvgGenerator
from qgis.PyQt.QtWidgets import (
    QDialog, QWidget, QVBoxLayout, QHBoxLayout, QTableWidget, QTabWidget,
    QTableWidgetItem, QPushButton, QComboBox, QLabel, QFileDialog,
    QMessageBox, QCheckBox, QDoubleSpinBox, QGroupBox, QGridLayout,
    QFrame
)

from qgis.core import (
    QgsProject, QgsRasterLayer, QgsVectorLayer, QgsFeature, QgsField,
    QgsGeometry, QgsPointXY, QgsWkbTypes, QgsRectangle,
    QgsCoordinateTransform, QgsCategorizedSymbolRenderer, QgsRendererCategory,
    QgsLineSymbol, QgsPalLayerSettings, QgsTextFormat, QgsVectorLayerSimpleLabeling,
    QgsUnitTypes, QgsSingleSymbolRenderer, QgsGraduatedSymbolRenderer,
    QgsVectorFileWriter,
    QgsMapLayer  # <-- AGGIUNTO QUI PER IL FIX DELLO STILE
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
# CLASSE STEREOPLOT DIALOG (INTEGRATA CON DISEGNO GRAFICO)
# =====================================================================
class StereoplotCanvas(QWidget):
    """
    Widget personalizzato avanzato per lo stereoplot (Proiezione di Schmidt Equiarea).
    Supporta Poli, Piani, Density Contour (Kamb/Kalsberg style) e Beta Axis.
    """
    def __init__(self, bedding_data, parent=None):
        super().__init__(parent)
        self.bedding_data = bedding_data  # Lista di tuple: (feature_id, azimuth, dip)
        self.setMinimumSize(400, 400)
        
        # Stati di visualizzazione controllati dai pulsanti
        self.show_poles = True
        self.show_planes = False
        self.show_contour = False
        self.show_beta = False
        
        self.beta_axis = None  # Memorizza (trend, plunge)
        self._contour_image = None
        
        # Precalcola i vettori 3D dei poli (nella semisfera inferior: z <= 0)
        self.pole_vectors = []
        for fid, az, dip in self.bedding_data:
            dip_dir_rad = math.radians(az)
            dip_rad = math.radians(dip)
            x = -math.sin(dip_rad) * math.sin(dip_dir_rad)
            y = -math.sin(dip_rad) * math.cos(dip_dir_rad)
            z = -math.cos(dip_rad)
            self.pole_vectors.append((x, y, z))

    def update_contour_grid(self, r):
        """Genera una heatmap di densità basata su un cerchio di conteggio ad area costante (5%)."""
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
        """Calcola l'asse Beta (Fold Axis) globale tramite l'analisi tensoriale."""
        if len(self.pole_vectors) < 2:
            self.beta_axis = None
            return

        m11, m12, m13 = 0.0, 0.0, 0.0
        m22, m23, m33 = 0.0, 0.0, 0.0

        for x, y, z in self.pole_vectors:
            m11 += x * x
            m12 += x * y
            m13 += x * z
            m22 += y * y
            m23 += y * z
            m33 += z * z

        A = [
            [m11, m12, m13],
            [m12, m22, m23],
            [m13, m23, m33]
        ]

        V = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        max_iterations = 50
        for _ in range(max_iterations):
            row, col = 0, 1
            max_val = abs(A[0][1])
            if abs(A[0][2]) > max_val:
                row, col = 0, 2
                max_val = abs(A[0][2])
            if abs(A[1][2]) > max_val:
                row, col = 1, 2
                max_val = abs(A[1][2])

            if max_val < 1e-9:
                break

            phi = 0.5 * math.atan2(2.0 * A[row][col], A[row][row] - A[col][col])
            c = math.cos(phi)
            s = math.sin(phi)

            ar_r = c * c * A[row][row] + 2.0 * s * c * A[row][col] + s * s * A[col][col]
            ac_c = s * s * A[row][row] - 2.0 * s * c * A[row][col] + c * c * A[col][col]
            A[row][row] = ar_r
            A[col][col] = ac_c
            A[row][col] = 0.0
            A[col][row] = 0.0

            other = 3 - row - col
            a_ro = c * A[row][other] + s * A[col][other]
            a_co = -s * A[row][other] + c * A[col][other]
            A[row][other] = A[other][row] = a_ro
            A[col][other] = A[other][col] = a_co

            for k in range(3):
                v_kr = c * V[k][row] + s * V[k][col]
                v_kc = -s * V[k][row] + c * V[k][col]
                V[k][row] = v_kr
                V[k][col] = v_kc

        eigenvalues = [A[0][0], A[1][1], A[2][2]]
        min_idx = eigenvalues.index(min(eigenvalues))
        bx = V[0][min_idx]
        by = V[1][min_idx]
        bz = V[2][min_idx]

        if bz > 0:
            bx, by, bz = -bx, -by, -bz

        horiz_dist = math.hypot(bx, by)
        plunge = math.degrees(math.atan2(abs(bz), horiz_dist))
        trend = math.degrees(math.atan2(bx, by)) % 360.0

        self.beta_axis = (trend, plunge)

    def draw_stereonet(self, painter, width, height):
        """Metodo centrale di rendering condiviso tra paintEvent ed esportazione SVG."""
        cx = width / 2
        cy = height / 2
        r = min(cx, cy) - 25

        # --- DENSITY CONTOUR BACKGROUND ---
        if self.show_contour:
            if self._contour_image is None or self._contour_image.width() != int(r*2):
                self.update_contour_grid(r)
            if self._contour_image:
                painter.drawImage(QPointF(cx - r, cy - r), self._contour_image)

        pen_grid = QPen(QColor(160, 160, 160), 1, Qt.PenStyle.DashLine)
        pen_border = QPen(Qt.GlobalColor.black, 2)
        
        # 1. Cerchio Primitivo ed Assi
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
            if math.isclose(h_mag, 0.0):
                return cx, cy
            px = cx + r_proj * (vx / h_mag)
            py = cy - r_proj * (vy / h_mag)
            return px, py

        # 2. PLOT DEI PIANI (Grandi Cerchi)
        if self.show_planes:
            pen_plane = QPen(QColor(50, 50, 220, 140), 1.5)
            painter.setPen(pen_plane)
            for fid, az, dip in self.bedding_data:
                dip_dir_rad = math.radians(az)
                dip_rad = math.radians(dip)
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

        # 3. PLOT DEI POLI
        if self.show_poles:
            pen_pole = QPen(Qt.GlobalColor.red, 5)
            painter.setPen(pen_pole)
            for x, y, z in self.pole_vectors:
                px, py = project_vector(x, y, z)
                painter.drawEllipse(QPointF(px, py), 3, 3)

        # 4. PLOT BETA AXIS
        if self.show_beta and self.beta_axis:
            b_trend, b_plunge = self.beta_axis
            tr_rad = math.radians(b_trend)
            pl_rad = math.radians(b_plunge)
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
    """
    Finestra interattiva dello stereoplot con barra degli strumenti a 4 pulsanti dedicati
    e opzione per l'esportazione vettoriale in SVG.
    """
    def __init__(self, parent, bedding_data):
        super().__init__(parent)
        self.setWindowTitle("Structural Stereoplot (Schmidt Projection)")
        self.resize(500, 650)
        self.bedding_data = bedding_data  
        
        layout = QVBoxLayout(self)
        
        # --- BARRA DEI PULSANTI (STRUMENTI DI ANALISI) ---
        grid_buttons = QGridLayout()
        
        self.btn_poles = QPushButton("1. Plot Poles")
        self.btn_poles.setCheckable(True)
        self.btn_poles.setChecked(True)
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
        
        # Area Grafica del Canvas
        self.stereoplot_canvas = StereoplotCanvas(self.bedding_data, self)
        layout.addWidget(self.stereoplot_canvas, stretch=1)
        
        # Etichetta informativa dinamica per i calcoli strutturali
        self.lbl_info = QLabel(f"Dataset: {len(self.bedding_data)} elementi misurati.")
        self.lbl_info.setStyleSheet("font-weight: bold; font-size: 10pt; color: #004d40; margin-left: 5px;")
        layout.addWidget(self.lbl_info)
        
        # --- PULSANTI DI AZIONE IN BASSO ---
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
        self.stereoplot_canvas.show_poles = checked
        self.stereoplot_canvas.update()

    def _toggle_planes(self, checked):
        self.stereoplot_canvas.show_planes = checked
        self.stereoplot_canvas.update()

    def _toggle_contour(self, checked):
        self.stereoplot_canvas.show_contour = checked
        self.stereoplot_canvas._contour_image = None 
        self.stereoplot_canvas.update()

    def _toggle_beta(self, checked):
        self.stereoplot_canvas.show_beta = checked
        if checked:
            self.stereoplot_canvas.compute_beta_axis()
            res = self.stereoplot_canvas.beta_axis
            if res:
                self.lbl_info.setText(f"Beta Axis (Fold Axis): Trend {res[0]:.1f}° / Plunge {res[1]:.1f}°")
            else:
                self.lbl_info.setText("Insufficient data number.")
        else:
            self.lbl_info.setText(f"Dataset: {len(self.bedding_data)} measured items.")
        self.stereoplot_canvas.update()

    def _export_stereoplot_to_svg(self):
        """Esporta la grafica corrente dello stereoplot in formato SVG."""
        path, _ = QFileDialog.getSaveFileName(self, "Export as SVG", "stereoplot.svg", "Scalable Vector Graphics (*.svg)")
        if not path:
            return

        # Dimensione quadrata fissa standard per l'esportazione del diagramma
        export_size = QSize(600, 600)
        
        generator = QSvgGenerator()
        generator.setFileName(path)
        generator.setSize(export_size)
        generator.setViewBox(QRectF(0, 0, export_size.width(), export_size.height()))
        generator.setTitle("Structural Stereoplot - Schmidt Projection")
        generator.setDescription(f"Generated from QGIS Parallel Folds Tool. Items count: {len(self.bedding_data)}")

        painter = QPainter(generator)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        
        # Sfondo bianco per il canvas SVG
        painter.fillRect(QRectF(0, 0, export_size.width(), export_size.height()), Qt.GlobalColor.white)
        
        # Invocazione del metodo di rendering sul generatore SVG
        self.stereoplot_canvas.draw_stereonet(painter, export_size.width(), export_size.height())
        painter.end()
        
        QMessageBox.information(self, "Exported", f"Stereoplot saved in SVG:\n{path}")


class FoldWindow(QDialog):
    
    def _on_tab_changed(self, index):
        if index != 1:
            if self.btn_start.isChecked():
                self.btn_start.setChecked(False)
            if hasattr(self, 'btn_edit_saved') and self.btn_edit_saved.isChecked():
                self.btn_edit_saved.setChecked(False)
        
    def __init__(self, iface):
        super().__init__(iface.mainWindow())
        self.iface = iface
        
        self.linea_separatrice_profile = QFrame()
        self.linea_separatrice_profile.setFrameShape(QFrame.Shape.HLine)
        self.linea_separatrice_profile.setFrameShadow(QFrame.Shadow.Sunken)
        
        self.linea_separatrice_folding = QFrame()
        self.linea_separatrice_folding.setFrameShape(QFrame.Shape.HLine)
        self.linea_separatrice_folding.setFrameShadow(QFrame.Shadow.Sunken)
        
        self.setWindowTitle("Parallel Folds & Structural Tool")
        self.resize(1150, 850)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)

        # --- MEMORIA DI PERSISTENZA DEI LAYER ---
        self.dem_layer = None
        self.trace_layer = None
        self.projected_structural_dips = []

        self.tool = None
        self.layers_thickness = [20.0, 20.0, 20.0]
        self.layers_colors = [DEFAULT_COLORS[i % len(DEFAULT_COLORS)] for i in range(16)]
        
        self.digitized_index = 0  
        self._canvas_layers = []
        self._output_layer = None
        self._ghost_layer = None  
        
        # Layer di memoria indipendenti
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

        # Layout orizzontale principale
        window_layout = QHBoxLayout()
        self.setLayout(window_layout)

        # =====================================================================
        # 1. PANNELLO DI SINISTRA
        # =====================================================================
        self.left_panel = QWidget()
        self.left_panel.setFixedWidth(340)
        left_layout = QVBoxLayout(self.left_panel)
        left_layout.setContentsMargins(0, 0, 0, 0)
        window_layout.addWidget(self.left_panel)

        # Widget dei Tab
        self.main_tabs = QTabWidget()
        left_layout.addWidget(self.main_tabs, stretch=1)

        self.tab_profile_extraction = QWidget()
        self.tab_folding_digitization = QWidget()
        self.tab_import_export = QWidget() # Nuovo tab per Import / Export
        
        self.main_tabs.addTab(self.tab_profile_extraction, "1. Data Extraction & Profile")
        self.main_tabs.addTab(self.tab_folding_digitization, "2. Parallel Fold Digitization")
        self.main_tabs.addTab(self.tab_import_export, "3. Import / Export")

        # Configurazione dei singoli tab
        self._setup_tab_profile()
        self._setup_tab_folding()
        self._setup_tab_import_export()

        # -----------------------------------------------------------------
        # PANNELLO VISIBILITÀ
        # -----------------------------------------------------------------
        self.visibility_container = QWidget()
        visibility_layout = QVBoxLayout(self.visibility_container)
        visibility_layout.setContentsMargins(15, 10, 15, 15)
        visibility_layout.setSpacing(8)

        # Intestazione
        lbl_visibility = QLabel("Layers visibility")
        lbl_visibility.setStyleSheet("font-weight: bold; font-size: 11px; color: #333333;")
        visibility_layout.addWidget(lbl_visibility)

        # Linea separatrice sotto l'intestazione
        line = QFrame()
        line.setFrameShape(QFrame.Shape.HLine)
        line.setFrameShadow(QFrame.Shadow.Sunken)
        visibility_layout.addWidget(line)

        # Definizione dei 5 checkbox
        self.chk_profile = QCheckBox("XY axes")
        self.chk_strikes = QCheckBox("Projected bedding")
        self.chk_traces = QCheckBox("Projected traces")
        self.chk_intersections = QCheckBox("Traces intersection")
        self.chk_parallel_folds = QCheckBox("Parallel Fold")

        # Impostiamo lo stato iniziale a True
        self.chk_profile.setChecked(True)
        self.chk_strikes.setChecked(True)
        self.chk_traces.setChecked(True)
        self.chk_intersections.setChecked(True)
        self.chk_parallel_folds.setChecked(True)

        # Aggiunta dei checkbox al layout verticale
        visibility_layout.addWidget(self.chk_profile)
        visibility_layout.addWidget(self.chk_strikes)
        visibility_layout.addWidget(self.chk_traces)
        visibility_layout.addWidget(self.chk_intersections)
        visibility_layout.addWidget(self.chk_parallel_folds)

        # Connessione dei segnali alle funzioni di toggle
        self.chk_profile.toggled.connect(lambda checked: self.toggle_layer_visibility("profile", checked))
        self.chk_strikes.toggled.connect(lambda checked: self.toggle_layer_visibility("strikes", checked))
        self.chk_traces.toggled.connect(lambda checked: self.toggle_layer_visibility("traces", checked))
        self.chk_intersections.toggled.connect(lambda checked: self.toggle_layer_visibility("intersections", checked))
        self.chk_parallel_folds.toggled.connect(lambda checked: self.toggle_layer_visibility("parallel_folds", checked))

        # Aggiungiamo il contenitore della visibilità in fondo al pannello di sinistra
        left_layout.addWidget(self.visibility_container)

        # =====================================================================
        # 2. CANVAS DI QGIS
        # =====================================================================
        self.map_canvas = QgsMapCanvas()
        self.map_canvas.setCanvasColor(QColor(255, 255, 255))
        self.map_canvas.setDestinationCrs(QgsProject.instance().crs())
        self.map_canvas.setExtent(QgsRectangle(0, 0, 1000, 1000))
        self.map_canvas.enableAntiAliasing(True)
        window_layout.addWidget(self.map_canvas, stretch=1)

        # Connessione dei segnali del Progetto QGIS
        QgsProject.instance().layersAdded.connect(self._populate_layers_combos)
        QgsProject.instance().layersRemoved.connect(self._populate_layers_combos)
        
        # Connessione del cambio Tab
        self.main_tabs.currentChanged.connect(self._on_tab_changed)

        # Aggiornamento dell'interfaccia
        self._refresh_table()
        self._refresh_combo()

    def toggle_layer_visibility(self, layer_key, visible):
        layer = None
        if layer_key == "profile":
            layer = getattr(self, "_profile_layer", None)
        elif layer_key == "strikes":
            layer = getattr(self, "_strikes_layer", None) 
        elif layer_key == "traces":
            layer = getattr(self, "_traces_layer", None)   
        elif layer_key == "intersections":
            layer = getattr(self, "_intersections_layer", None) 
        elif layer_key == "parallel_folds":
            layer = getattr(self, "_output_layer", None)

        if layer is None:
            return

        # Aggiorna nel Layer Tree di QGIS
        root = QgsProject.instance().layerTreeRoot()
        layer_node = root.findLayer(layer.id())
        if layer_node:
            layer_node.setItemVisibilityChecked(visible)

        # Aggiorna il canvas interno
        self._update_canvas_visibility()
        
    def _update_canvas_visibility(self):
        if not hasattr(self, "map_canvas") or not hasattr(self, "_canvas_layers"):
            return

        visible_layers = []
        for layer in self._canvas_layers:
            if not layer:
                continue
                
            if layer == getattr(self, "_profile_layer", None):
                if not self.chk_profile.isChecked():
                    continue
            elif layer == getattr(self, "_strikes_layer", None):
                if not self.chk_strikes.isChecked():
                    continue
            elif layer == getattr(self, "_traces_layer", None):
                if not self.chk_traces.isChecked():
                    continue
            elif layer == getattr(self, "_intersections_layer", None):
                if not self.chk_intersections.isChecked():
                    continue
            elif layer == getattr(self, "_output_layer", None):
                if not self.chk_parallel_folds.isChecked():
                    continue

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

        # TAB 1: DEM E TRACCIA
        sub_layout_1 = QVBoxLayout(self.sub_tab_dem_trace)
        sub_layout_1.addWidget(QLabel("Select Digital Elevation Model (DEM):"))
        self.combo_dem = QComboBox()
        sub_layout_1.addWidget(self.combo_dem)

        sub_layout_1.addWidget(QLabel("Select Section Line (Vector Layer):"))
        self.combo_trace = QComboBox()
        sub_layout_1.addWidget(self.linea_separatrice_profile)
        sub_layout_1.addWidget(self.combo_trace)

        sub_layout_1.addWidget(QLabel("Sampling Interval (Points count):"))
        self.spin_sampling = QDoubleSpinBox()
        self.spin_sampling.setRange(10.0, 5000.0)
        self.spin_sampling.setValue(200.0)
        sub_layout_1.addWidget(self.spin_sampling)

        self.btn_generate_profile = QPushButton("Extract & Generate Profile")
        self.btn_generate_profile.setStyleSheet("background-color: #2e7d32; color: white; font-weight: bold; font-size: 11pt;")
        self.btn_generate_profile.clicked.connect(self._generate_topographic_profile)
        sub_layout_1.addWidget(self.btn_generate_profile)
        sub_layout_1.addStretch()

        # TAB 2: APPARENT DIP
        layout_2 = QVBoxLayout(self.sub_tab_2)
        layout_2.addWidget(QLabel("Select Structural Points Layer:"))
        self.combo_struct_layer = QComboBox()
        layout_2.addWidget(self.combo_struct_layer)
        self.combo_struct_layer.currentIndexChanged.connect(self._on_struct_layer_changed)

        layout_2.addWidget(QLabel("Layer Orientation Format:"))
        self.combo_orient_format = QComboBox()
        self.combo_orient_format.addItems(["Dip Direction / Dip", "Strike / Dip (Right-Hand Rule)"])
        layout_2.addWidget(self.combo_orient_format)

        layout_dir_dip = QHBoxLayout()
        vbox_dir = QVBoxLayout()
        self.lbl_dir_field = QLabel("Dip Direction Field:")
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

        group_proj = QGroupBox("Custom Projection Vector (Calculated automatically if left 90/0)")
        layout_proj = QGridLayout()
        layout_proj.addWidget(QLabel("Projection Trend (0-360°):"), 0, 0)
        self.spin_proj_trend = QDoubleSpinBox()
        self.spin_proj_trend.setRange(0.0, 360.0)
        self.spin_proj_trend.setValue(90.0)
        layout_proj.addWidget(self.spin_proj_trend, 0, 1)
        layout_proj.addWidget(QLabel("Projection Plunge (0-90°):"), 1, 0)
        self.spin_proj_plunge = QDoubleSpinBox()
        self.spin_proj_plunge.setRange(0.0, 90.0)
        self.spin_proj_plunge.setValue(0.0)
        layout_proj.addWidget(self.spin_proj_plunge, 1, 1)
        group_proj.setLayout(layout_proj)
        layout_2.addWidget(group_proj)

        layout_h_params = QHBoxLayout()
        vbox_buf = QVBoxLayout()
        vbox_buf.addWidget(QLabel("Max Search Distance (m):"))
        self.spin_dip_buffer = QDoubleSpinBox()
        self.spin_dip_buffer.setRange(1.0, 100000.0)
        self.spin_dip_buffer.setValue(1000.0)
        vbox_buf.addWidget(self.spin_dip_buffer)
        layout_h_params.addLayout(vbox_buf)

        vbox_len = QVBoxLayout()
        vbox_len.addWidget(QLabel("Symbol Length (m):"))
        self.spin_dip_symbol_len = QDoubleSpinBox()
        self.spin_dip_symbol_len.setRange(5.0, 2000.0)
        self.spin_dip_symbol_len.setValue(50.0)
        vbox_len.addWidget(self.spin_dip_symbol_len)
        layout_h_params.addLayout(vbox_len)
        layout_2.addLayout(layout_h_params)

        self.btn_project_dips = QPushButton("Project Apparent Dips")
        self.btn_project_dips.setStyleSheet("background-color: #ef6c00; color: white; font-weight: bold;")
        self.btn_project_dips.clicked.connect(self._compute_apparent_dips)
        layout_2.addWidget(self.btn_project_dips)
        
        # --- PULSANTE STEREOPLOT INSERITO QUI ---
        self.btn_open_stereoplot = QPushButton("Open Stereoplot Window")
        self.btn_open_stereoplot.setStyleSheet("background-color: #00897b; color: white; font-weight: bold;")
        self.btn_open_stereoplot.clicked.connect(self._open_stereoplot)
        layout_2.addWidget(self.btn_open_stereoplot)
        
        layout_2.addStretch()

        # TAB 3: TRACE INTERSECTION
        layout_project_traces = QVBoxLayout(self.sub_tab_project_traces)
        layout_project_traces.addWidget(QLabel("Project Traces onto Profile:"))
        
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

        # TAB 4: PROJECT TRACES
        layout_project_traces_four = QVBoxLayout(self.sub_tab_project_traces_four)
        
        layout_project_traces_four.addWidget(QLabel("Select Traces Layer (Lines/Polygons):"))
        self.combo_pt4_traces_layer = QComboBox()
        layout_project_traces_four.addWidget(self.combo_pt4_traces_layer)
        self.combo_pt4_traces_layer.currentIndexChanged.connect(self._on_pt4_traces_layer_changed)

        layout_project_traces_four.addWidget(QLabel("Select Label Attribute Field (Optional):"))
        self.combo_pt4_label_field = QComboBox()
        layout_project_traces_four.addWidget(self.combo_pt4_label_field)

        self.chk_pt4_selected_only = QCheckBox("Project selected features only in QGIS")
        layout_project_traces_four.addWidget(self.chk_pt4_selected_only)

        group_proj_four = QGroupBox("Projection Vector")
        layout_proj_four = QGridLayout()
        layout_proj_four.addWidget(QLabel("Projection Trend (0-360°):"), 0, 0)
        self.spin_pt4_trend = QDoubleSpinBox()
        self.spin_pt4_trend.setRange(0.0, 360.0)
        self.spin_pt4_trend.setValue(90.0)
        layout_proj_four.addWidget(self.spin_pt4_trend, 0, 1)
        
        layout_proj_four.addWidget(QLabel("Projection Plunge (0-90°):"), 1, 0)
        self.spin_pt4_plunge = QDoubleSpinBox()
        self.spin_pt4_plunge.setRange(0.0, 90.0)
        self.spin_pt4_plunge.setValue(0.0)
        layout_proj_four.addWidget(self.spin_pt4_plunge, 1, 1)
        group_proj_four.setLayout(layout_proj_four)
        layout_project_traces_four.addWidget(group_proj_four)

        layout_pt4_dist = QHBoxLayout()
        layout_pt4_dist.addWidget(QLabel("Max Search Distance (m):"))
        self.spin_pt4_buffer = QDoubleSpinBox()
        self.spin_pt4_buffer.setRange(1.0, 100000.0)
        self.spin_pt4_buffer.setValue(1000.0)
        layout_pt4_dist.addWidget(self.spin_pt4_buffer)
        layout_project_traces_four.addLayout(layout_pt4_dist)

        self.btn_pt4_project = QPushButton("Project Traces Nodes")
        self.btn_pt4_project.setStyleSheet("background-color: #6a1b9a; color: white; font-weight: bold;")
        self.btn_pt4_project.clicked.connect(self._compute_projected_traces_nodes)
        layout_project_traces_four.addWidget(self.btn_pt4_project)

        self.btn_pt4_select_single = QPushButton("Select & Delete Single Projected Line")
        self.btn_pt4_select_single.setCheckable(True)
        self.btn_pt4_select_single.setStyleSheet("background-color: #0288d1; color: white; font-weight: bold;")
        self.btn_pt4_select_single.clicked.connect(self._toggle_select_proj_tool)
        layout_project_traces_four.addWidget(self.btn_pt4_select_single)

        self.btn_pt4_delete_selected = QPushButton("Delete Selected Projected Line")
        self.btn_pt4_delete_selected.setStyleSheet("background-color: #e53935; color: white; font-weight: bold;")
        self.btn_pt4_delete_selected.clicked.connect(self._delete_selected_projected_line)
        layout_project_traces_four.addWidget(self.btn_pt4_delete_selected)

        self.btn_pt4_clear = QPushButton("Clear All Projected Traces")
        self.btn_pt4_clear.setStyleSheet("background-color: #c62828; color: white; font-weight: bold;")
        self.btn_pt4_clear.clicked.connect(self._clear_projected_traces_nodes)
        layout_project_traces_four.addWidget(self.btn_pt4_clear)

        layout_project_traces_four.addStretch()
        self._populate_layers_combos()

    # =====================================================================
    # METODO: APERTURA DIALOG E GESTIONE DATI STEREOPLOT
    # =====================================================================
    def _open_stereoplot(self):
        """
        Filtra i punti di bedding compresi nella distanza di buffer e nella estensione 
        della linea di sezione, e apre la nuova finestra 'stereoplot'.
        """
        if not self.trace_layer:
            QMessageBox.warning(self, "Missing trace", "Extract the topographic profile first to define the section line.")
            return

        struct_layer_id = self.combo_struct_layer.currentData()
        struct_layer = QgsProject.instance().mapLayer(struct_layer_id) if struct_layer_id else None
        dir_field = self.combo_dir_field.currentData()
        dip_field = self.combo_dip_field.currentData()

        if not struct_layer or not dir_field or not dip_field:
            QMessageBox.warning(self, "Missing Selection", "Verify that the structural layer and Azimuth/Dip fields are correctly selected.")
            return

        section_features = list(self.trace_layer.getFeatures())
        if not section_features:
            return

        section_geom = QgsGeometry(section_features[0].geometry())
        section_geom.get().dropZValue()
        section_geom.get().dropMValue()

        # Coordinate transform per lavorare nello stesso sistema del DEM/Sezione
        xform_section_to_dem = QgsCoordinateTransform(self.trace_layer.crs(), self.dem_layer.crs(), QgsProject.instance())
        section_geom_dem_crs = QgsGeometry(section_geom)
        section_geom_dem_crs.transform(xform_section_to_dem)

        if section_geom_dem_crs.isMultipart():
            parts = section_geom_dem_crs.asMultiPolyline()
            polyline = parts[0] if parts else []
        else:
            polyline = section_geom_dem_crs.asPolyline()
            
        if len(polyline) < 2: 
            return

        buffer_dist = self.spin_dip_buffer.value()
        format_idx = self.combo_orient_format.currentIndex()
        section_length = section_geom_dem_crs.length()

        # Parametri di proiezione vettoriale vettorializzati per il calcolo dell'intersezione reale
        proj_trend = self.spin_proj_trend.value()
        proj_plunge = self.spin_proj_plunge.value()
        alpha_proj = math.radians(proj_trend)
        beta_proj = math.radians(proj_plunge)

        px = math.sin(alpha_proj) * math.cos(beta_proj)
        py = math.cos(alpha_proj) * math.cos(beta_proj)

        p_start = polyline[0]
        p_end = polyline[-1]
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
            if geom.isEmpty():
                continue

            pt_dem_geom = QgsGeometry(geom)
            pt_dem_geom.transform(xform_struct_to_dem)
            pt_dem = pt_dem_geom.asPoint()

            # 1. Filtro preliminare Bounding Box e Distanza 2D massima dalla traccia
            if not section_bbox.contains(pt_dem): 
                continue
            if section_geom_dem_crs.distance(pt_dem_geom) > buffer_dist: 
                continue

            # 2. Calcolo del punto proiettato sulla sezione lungo il vettore di structural trend
            if math.isclose(denom, 0.0):
                continue

            t = -(A * pt_dem.x() + B * pt_dem.y() + C) / denom
            x_inter = pt_dem.x() + t * px
            y_inter = pt_dem.y() + t * py

            proj_point = QgsPointXY(x_inter, y_inter)
            proj_geom_pt = QgsGeometry.fromPointXY(proj_point)

            # 3. Verifica se la proiezione cade longitudinalmente dentro l'estensione fisica della sezione
            dist_local = section_geom_dem_crs.lineLocatePoint(proj_geom_pt)
            if dist_local < 0.0 or dist_local > section_length:
                continue

            # Se supera tutti i filtri di prossimità, estrae e normalizza l'orientazione
            try:
                raw_dir = float(feature[dir_field])
                dip = float(feature[dip_field])
                
                azimuth = raw_dir % 360.0 if format_idx == 0 else (raw_dir + 90.0) % 360.0
                bedding_within_buffer.append((feature.id(), azimuth, dip))
            except (ValueError, TypeError):
                continue

        if not bedding_within_buffer:
            QMessageBox.information(self, "No Data", f"No bedding found within the limits.")
            return

        # Apre la finestra dello stereoplot passando solo i dati validati e coerenti
        self.stereoplot_window = StereoplotDialog(self, bedding_within_buffer)
        self.stereoplot_window.show()
        

    def _populate_layers_combos(self):
        if getattr(self, '_block_combo_updates', False):
            return

        self.combo_dem.blockSignals(True)
        self.combo_trace.blockSignals(True)
        self.combo_intersection_traces.blockSignals(True)
        if hasattr(self, 'combo_pt4_traces_layer'):
            self.combo_pt4_traces_layer.blockSignals(True)
        if hasattr(self, 'combo_struct_layer'):
            self.combo_struct_layer.blockSignals(True)

        self.combo_dem.clear()
        self.combo_trace.clear()
        self.combo_intersection_traces.clear()
        if hasattr(self, 'combo_pt4_traces_layer'):
            self.combo_pt4_traces_layer.clear()
        if hasattr(self, 'combo_struct_layer'):
            self.combo_struct_layer.clear()

        excluded_keywords = {
            "topographic profile with grid",
            "projected bedding",
            "projected traces",
            "traces intersection",
            "profile base elements (ghost)",
            "parallel & structural layers"
        }

        for layer_id, layer in QgsProject.instance().mapLayers().items():
            layer_name_lower = layer.name().lower()
            
            if any(kw in layer_name_lower for kw in excluded_keywords):
                continue

            if isinstance(layer, QgsRasterLayer):
                if "background" in layer_name_lower or "calibrated" in layer_name_lower:
                    continue
                self.combo_dem.addItem(layer.name(), layer_id)
            elif isinstance(layer, QgsVectorLayer):
                g_type = layer.geometryType()
                if g_type == 0:
                    if hasattr(self, 'combo_struct_layer'):
                        self.combo_struct_layer.addItem(layer.name(), layer_id)
                elif g_type == 1:
                    self.combo_trace.addItem(layer.name(), layer_id)
                    self.combo_intersection_traces.addItem(layer.name(), layer_id)
                    if hasattr(self, 'combo_pt4_traces_layer'):
                        self.combo_pt4_traces_layer.addItem(layer.name(), layer_id)
                elif g_type == 2:
                    self.combo_intersection_traces.addItem(layer.name(), layer_id)
                    if hasattr(self, 'combo_pt4_traces_layer'):
                        self.combo_pt4_traces_layer.addItem(layer.name(), layer_id)

        self.combo_dem.blockSignals(False)
        self.combo_trace.blockSignals(False)
        self.combo_intersection_traces.blockSignals(False)
        if hasattr(self, 'combo_pt4_traces_layer'):
            self.combo_pt4_traces_layer.blockSignals(False)
        if hasattr(self, 'combo_struct_layer'):
            self.combo_struct_layer.blockSignals(False)
        self._on_traces_layer_changed()
        self._on_pt4_traces_layer_changed()
        self._on_struct_layer_changed()
        

    def _on_traces_layer_changed(self):
        self.combo_trace_label_field.clear()
        self.combo_trace_label_field.addItem("None", None)
        l_id = self.combo_intersection_traces.currentData()
        layer = QgsProject.instance().mapLayer(l_id) if l_id else None
        if layer and isinstance(layer, QgsVectorLayer):
            for f in layer.fields():
                self.combo_trace_label_field.addItem(f.name(), f.name())

    def _on_pt4_traces_layer_changed(self):
        if not hasattr(self, 'combo_pt4_label_field'): return
        self.combo_pt4_label_field.clear()
        self.combo_pt4_label_field.addItem("None (No Label)", None)
        l_id = self.combo_pt4_traces_layer.currentData()
        layer = QgsProject.instance().mapLayer(l_id) if l_id else None
        if layer and isinstance(layer, QgsVectorLayer):
            for f in layer.fields():
                self.combo_pt4_label_field.addItem(f.name(), f.name())

    def _on_struct_layer_changed(self):
        if not hasattr(self, 'combo_struct_layer'): return
        self.combo_dir_field.clear()
        self.combo_dip_field.clear()
        l_id = self.combo_struct_layer.currentData()
        layer = QgsProject.instance().mapLayer(l_id) if l_id else None
        if layer and isinstance(layer, QgsVectorLayer):
            for f in layer.fields():
                if f.isNumeric():
                    self.combo_dir_field.addItem(f.name(), f.name())
                    self.combo_dip_field.addItem(f.name(), f.name())

    def _generate_topographic_profile(self):
        dem_layer_id = self.combo_dem.currentData()
        trace_layer_id = self.combo_trace.currentData()

        if not dem_layer_id or not trace_layer_id:
            QMessageBox.warning(self, "Missing Layers", "Please select both a DEM and a Trace layer.")
            return

        self.dem_layer = QgsProject.instance().mapLayer(dem_layer_id)
        self.trace_layer = QgsProject.instance().mapLayer(trace_layer_id)

        if not self.dem_layer or not self.trace_layer:
            QMessageBox.critical(self, "Error", "Selected layers are invalid.")
            return

        features = list(self.trace_layer.getFeatures())
        if not features:
            QMessageBox.critical(self, "Trace Empty", "The selected trace layer contains no lines.")
            return

        geom = QgsGeometry(features[0].geometry())
        geom.get().dropZValue()
        geom.get().dropMValue()

        xform = QgsCoordinateTransform(self.trace_layer.crs(), self.dem_layer.crs(), QgsProject.instance())
        geom_dem_crs = QgsGeometry(geom)
        geom_dem_crs.transform(xform)

        self.profile_length = geom_dem_crs.length()
        if self.profile_length <= 0:
            QMessageBox.critical(self, "Error", "The section line length is zero.")
            return

        num_samples = int(self.spin_sampling.value())
        self.topo_points = []
        dp = self.dem_layer.dataProvider()

        self.profile_z_min = 999999.0
        self.profile_z_max = -999999.0

        for i in range(num_samples):
            fraction = i / (num_samples - 1)
            dist_local = fraction * self.profile_length
            pt_geom = geom_dem_crs.interpolate(dist_local)
            if pt_geom.isEmpty(): continue
            pt_dem = pt_geom.asPoint()

            val, ok = dp.sample(pt_dem, 1)
            qz = val if (ok and not math.isnan(val)) else 0.0

            if qz < self.profile_z_min: self.profile_z_min = qz
            if qz > self.profile_z_max: self.profile_z_max = qz
            self.topo_points.append(QgsPointXY(dist_local, qz))

        if self._profile_layer:
            self._remove_canvas_layer(self._profile_layer)
            self._profile_layer = None

        crs_canvas = self.map_canvas.mapSettings().destinationCrs()
        self._profile_layer = QgsVectorLayer(f"LineString?crs={crs_canvas.authid()}", "Topographic Profile with Grid", "memory")
        provider = self._profile_layer.dataProvider()
        
        provider.addAttributes([
            QgsField("id", QMetaType.Type.Int), 
            QgsField("label", QMetaType.Type.QString),
            QgsField("color", QMetaType.Type.QString)
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
        
        text_format = QgsTextFormat()
        text_format.setFont(QFont("Arial"))
        text_format.setSizeUnit(QgsUnitTypes.RenderUnit.RenderPixels) 
        text_format.setSize(12.0)  
        text_format.setColor(QColor("#333333"))
        
        label_settings.setFormat(text_format)
        label_settings.placement = QgsPalLayerSettings.Placement.Line
        label_settings.placementFlags = QgsPalLayerSettings.LinePlacementFlags.BelowLine
        label_settings.obstacle = False                             
        label_settings.allowDegradedPlacement = True                
        label_settings.priority = 10                                
        label_settings.dist = 2.0 
        
        self._profile_layer.setLabeling(QgsVectorLayerSimpleLabeling(label_settings))
        self._profile_layer.setLabelsEnabled(True)

        self._add_canvas_layer(self._profile_layer, zoom_to_it=True)
        self._regenerate_ghost_layer()
        self._update_canvas_visibility()
        
    def _update_grid_and_axes(self):
        if not self._profile_layer: return
        provider = self._profile_layer.dataProvider()
        
        old_ids = [f.id() for f in self._profile_layer.getFeatures() if f["id"] in (2, 3)]
        if old_ids: 
            provider.deleteFeatures(old_ids)

        margin_y = (self.profile_z_max - self.profile_z_min) * 0.25 if (self.profile_z_max - self.profile_z_min) > 0 else 10.0
        axis_y_bottom = self.profile_z_min - margin_y
        axis_y_top = self.profile_z_max + margin_y

        len_x = self.profile_length
        len_y = axis_y_top - axis_y_bottom

        tick_size_comune = len_x / 50.0

        features_axes = []
        
        feat_ax = QgsFeature(self._profile_layer.fields())
        feat_ax.setGeometry(QgsGeometry.fromPolylineXY([QgsPointXY(0.0, axis_y_bottom), QgsPointXY(len_x, axis_y_bottom)]))
        feat_ax.setAttributes([2, "", ""])
        features_axes.append(feat_ax)

        feat_ay = QgsFeature(self._profile_layer.fields())
        feat_ay.setGeometry(QgsGeometry.fromPolylineXY([QgsPointXY(0.0, axis_y_bottom), QgsPointXY(0.0, axis_y_top)]))
        feat_ay.setAttributes([2, "", ""])
        features_axes.append(feat_ay)

        step_x = self._calculate_nice_step(len_x)
        tick_length_y = tick_size_comune
        
        curr_x = 0.0
        while curr_x <= len_x:
            feat_tick = QgsFeature(self._profile_layer.fields())
            feat_tick.setGeometry(QgsGeometry.fromPolylineXY([
                QgsPointXY(curr_x, axis_y_bottom), 
                QgsPointXY(curr_x, axis_y_bottom - tick_length_y)
            ]))
            feat_tick.setAttributes([3, f"{curr_x:g}", ""])
            features_axes.append(feat_tick)
            curr_x += step_x

        step_y = self._calculate_nice_step(len_y) / 0.5
        tick_length_x = tick_size_comune
        
        curr_y = math.ceil(axis_y_bottom / step_y) * step_y
        while curr_y <= axis_y_top:
            feat_tick_y = QgsFeature(self._profile_layer.fields())
            feat_tick_y.setGeometry(QgsGeometry.fromPolylineXY([
                QgsPointXY(0.0, curr_y), 
                QgsPointXY(-tick_length_x, curr_y)
            ]))
            feat_tick_y.setAttributes([3, f"{curr_y:g}", ""])
            features_axes.append(feat_tick_y)
            curr_y += step_y

        provider.addFeatures(features_axes)
        
        self._profile_layer.updateExtents()
        self._profile_layer.updateFields()
        self._profile_layer.triggerRepaint()
        
    def _compute_intersections(self):
        if not self.dem_layer or not self.trace_layer or not self._profile_layer:
            QMessageBox.warning(self, "No Profile", "Please extract the topographic profile first to initialize data memory.")
            return

        traces_layer_id = self.combo_intersection_traces.currentData()
        traces_layer = QgsProject.instance().mapLayer(traces_layer_id) if traces_layer_id else None
        if not traces_layer:
            QMessageBox.warning(self, "Missing Layer", "Please select a valid vector Traces layer.")
            return

        section_features = list(self.trace_layer.getFeatures())
        if not section_features: return

        section_geom = QgsGeometry(section_features[0].geometry())
        section_geom.get().dropZValue()
        section_geom.get().dropMValue()

        xform_section_to_dem = QgsCoordinateTransform(self.trace_layer.crs(), self.dem_layer.crs(), QgsProject.instance())
        section_geom_dem_crs = QgsGeometry(section_geom)
        section_geom_dem_crs.transform(xform_section_to_dem)

        xform_section_to_traces = QgsCoordinateTransform(self.trace_layer.crs(), traces_layer.crs(), QgsProject.instance())
        section_geom_traces_crs = QgsGeometry(section_geom)
        section_geom_traces_crs.transform(xform_section_to_traces)

        xform_traces_to_dem = QgsCoordinateTransform(traces_layer.crs(), self.dem_layer.crs(), QgsProject.instance())
        label_field = self.combo_trace_label_field.currentData()
        dp = self.dem_layer.dataProvider()

        intersections_found = []
        for feature in traces_layer.getFeatures():
            t_geom = feature.geometry()
            if t_geom.isEmpty(): continue
            t_geom.get().dropZValue()
            t_geom.get().dropMValue()

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

        if not intersections_found:
            QMessageBox.information(self, "No Results", "No intersections found between Traces and Profile.")
            return

        if self._intersections_layer:
            self._remove_canvas_layer(self._intersections_layer)
            self._intersections_layer = None

        crs_canvas = self.map_canvas.mapSettings().destinationCrs()
        self._intersections_layer = QgsVectorLayer(f"LineString?crs={crs_canvas.authid()}", "Traces Intersection", "memory")
        provider = self._intersections_layer.dataProvider()
        
        provider.addAttributes([
            QgsField("id", QMetaType.Type.Int), 
            QgsField("label", QMetaType.Type.QString),
            QgsField("color", QMetaType.Type.QString)
        ])
        self._intersections_layer.updateFields()

        margin_y = (self.profile_z_max - self.profile_z_min) * 0.05 if (self.profile_z_max - self.profile_z_min) > 0 else 10.0
        axis_y_bottom = self.profile_z_min - margin_y

        features_to_add = []
        for dist, qz, lbl in intersections_found:
            feat = QgsFeature(self._intersections_layer.fields())
            feat.setGeometry(QgsGeometry.fromPolylineXY([QgsPointXY(dist, axis_y_bottom), QgsPointXY(dist, qz)]))
            feat.setAttributes([4, lbl, ""])
            features_to_add.append(feat)

        provider.addFeatures(features_to_add)
        
        sym = QgsLineSymbol.createSimple({'color': '#d32f2f', 'width': '1.0', 'line_style': 'dash'})
        self._intersections_layer.setRenderer(QgsSingleSymbolRenderer(sym))

        label_settings = QgsPalLayerSettings()
        label_settings.fieldName = "label"
        text_format = QgsTextFormat()
        text_format.setFont(QFont("Arial"))
        text_format.setSizeUnit(QgsUnitTypes.RenderUnit.RenderPixels) 
        text_format.setSize(12.0)  
        text_format.setColor(QColor("#d32f2f"))
        label_settings.setFormat(text_format)
        label_settings.placement = QgsPalLayerSettings.Placement.Line
        label_settings.placementFlags = QgsPalLayerSettings.LinePlacementFlags.AboveLine
        label_settings.obstacle = False                             
        label_settings.allowDegradedPlacement = True                
        label_settings.priority = 8                                
        label_settings.dist = 2.0 
        self._intersections_layer.setLabeling(QgsVectorLayerSimpleLabeling(label_settings))
        self._intersections_layer.setLabelsEnabled(True)

        self._intersections_layer.updateExtents()
        self._add_canvas_layer(self._intersections_layer)
        self._update_canvas_visibility()
        
        QMessageBox.information(self, "Success", f"Projected {len(intersections_found)} intersection lines.")

    def _compute_apparent_dips(self):
        if not self.dem_layer or not self.trace_layer or not self._profile_layer:
            if hasattr(self, 'iface') and self.iface:
                self.iface.messageBar().pushMessage("Error", "Please extract the topographic profile first.", level=1, duration=3)
            return

        struct_layer_id = self.combo_struct_layer.currentData()
        struct_layer = QgsProject.instance().mapLayer(struct_layer_id) if struct_layer_id else None
        dir_field = self.combo_dir_field.currentData()
        dip_field = self.combo_dip_field.currentData()

        if not struct_layer or not dir_field or not dip_field:
            if hasattr(self, 'iface') and self.iface:
                self.iface.messageBar().pushMessage("Error", "Verify structural layer and numeric fields selections.", level=1, duration=3)
            return

        section_features = list(self.trace_layer.getFeatures())
        if not section_features: return

        section_geom = QgsGeometry(section_features[0].geometry())
        section_geom.get().dropZValue()
        section_geom.get().dropMValue()

        xform_section_to_dem = QgsCoordinateTransform(self.trace_layer.crs(), self.dem_layer.crs(), QgsProject.instance())
        section_geom_dem_crs = QgsGeometry(section_geom)
        section_geom_dem_crs.transform(xform_section_to_dem)

        if section_geom_dem_crs.isMultipart():
            parts = section_geom_dem_crs.asMultiPolyline()
            if parts:
                polyline = parts[0]
            else:
                polyline = []
        else:
            polyline = section_geom_dem_crs.asPolyline()
            
        if len(polyline) < 2: return

        dx_s = polyline[-1].x() - polyline[0].x()
        dy_s = polyline[-1].y() - polyline[0].y()
        section_azimuth = math.degrees(math.atan2(dx_s, dy_s)) % 360.0

        buffer_dist = self.spin_dip_buffer.value()
        symbol_len = self.spin_dip_symbol_len.value()
        format_idx = self.combo_orient_format.currentIndex()

        proj_trend = self.spin_proj_trend.value()
        proj_plunge = self.spin_proj_plunge.value()

        xform_struct_to_dem = QgsCoordinateTransform(struct_layer.crs(), self.dem_layer.crs(), QgsProject.instance())
        dp = self.dem_layer.dataProvider()
        
        self.projected_structural_dips = []
        section_length = section_geom_dem_crs.length()

        section_bbox = section_geom_dem_crs.boundingBox()
        section_bbox.grow(buffer_dist + 10.0)

        for feature in struct_layer.getFeatures():
            geom = feature.geometry()
            if geom.isEmpty(): continue

            pt_dem_geom = QgsGeometry(geom)
            pt_dem_geom.transform(xform_struct_to_dem)
            pt_dem = pt_dem_geom.asPoint()

            if not section_bbox.contains(pt_dem): continue
            if section_geom_dem_crs.distance(pt_dem_geom) > buffer_dist: continue

            try:
                raw_dir = float(feature[dir_field])
                real_dip = float(feature[dip_field])
            except (ValueError, TypeError):
                continue

            dip_dir = raw_dir % 360.0 if format_idx == 0 else (raw_dir + 90.0) % 360.0
            
            angle_diff = math.radians(dip_dir - section_azimuth)
            apparent_dip_rad = math.atan(math.tan(math.radians(real_dip)) * math.cos(angle_diff))
            apparent_dip_deg = math.degrees(apparent_dip_rad)

            alpha_proj = math.radians(proj_trend)
            beta_proj = math.radians(proj_plunge)

            px = math.sin(alpha_proj) * math.cos(beta_proj)
            py = math.cos(alpha_proj) * math.cos(beta_proj)
            pz = -math.sin(beta_proj)

            val_init, ok_init = dp.sample(pt_dem, 1)
            z_init = val_init if (ok_init and not math.isnan(val_init)) else 0.0

            p_start = polyline[0]
            p_end = polyline[-1]
            
            A = p_end.y() - p_start.y()
            B = p_start.x() - p_end.x()
            C = p_end.x() * p_start.y() - p_start.x() * p_end.y()

            denom = A * px + B * py
            if math.isclose(denom, 0.0): 
                continue  

            t = -(A * pt_dem.x() + B * pt_dem.y() + C) / denom

            x_inter = pt_dem.x() + t * px
            y_inter = pt_dem.y() + t * py
            qz_projected = z_init + t * pz  

            proj_point = QgsPointXY(x_inter, y_inter)
            proj_geom_pt = QgsGeometry.fromPointXY(proj_point)
            
            if not section_geom_dem_crs.intersects(proj_geom_pt.buffer(0.001, 2)):
                continue 

            dist_local = section_geom_dem_crs.lineLocatePoint(proj_geom_pt)
            if dist_local < 0.0 or dist_local > section_length:
                continue

            self.projected_structural_dips.append((dist_local, qz_projected, apparent_dip_deg))

        if not self.projected_structural_dips:
            if hasattr(self, 'iface') and self.iface:
                self.iface.messageBar().pushMessage("Warning", "No structural points fall inside the section line boundaries after vector projection.", level=1, duration=3)
            return

        if self._strikes_layer:
            self._remove_canvas_layer(self._strikes_layer)
            self._strikes_layer = None

        crs_canvas = self.map_canvas.mapSettings().destinationCrs()
        self._strikes_layer = QgsVectorLayer(f"LineString?crs={crs_canvas.authid()}", "Projected Bedding", "memory")
        provider = self._strikes_layer.dataProvider()
        provider.addAttributes([
            QgsField("id", QMetaType.Type.Int), 
            QgsField("label", QMetaType.Type.QString),
            QgsField("color", QMetaType.Type.QString)
        ])
        self._strikes_layer.updateFields()

        features_to_add = []
        half_len = symbol_len / 2.0

        for dist, qz_projected, app_dip in self.projected_structural_dips:
            angle_rad = math.radians(app_dip)
            dx_p = half_len * math.cos(angle_rad)
            dy_p = -half_len * math.sin(angle_rad)

            feat_line = QgsFeature(self._strikes_layer.fields())
            feat_line.setGeometry(QgsGeometry.fromPolylineXY([
                QgsPointXY(dist - dx_p, qz_projected - dy_p), 
                QgsPointXY(dist + dx_p, qz_projected + dy_p)
            ]))
            feat_line.setAttributes([5, f"{abs(app_dip):.1f}°", ""])
            features_to_add.append(feat_line)

            feat_dot = QgsFeature(self._strikes_layer.fields())
            feat_dot.setGeometry(QgsGeometry.fromPolylineXY([
                QgsPointXY(dist - symbol_len * 0.02, qz_projected), 
                QgsPointXY(dist + symbol_len * 0.02, qz_projected)
            ]))
            feat_dot.setAttributes([6, "", ""])
            features_to_add.append(feat_dot)

        provider.addFeatures(features_to_add)

        categories = []
        categories.append(QgsRendererCategory(5, QgsLineSymbol.createSimple({'color': '#2e7d32', 'width': '1.5'}), "Apparent Dip"))
        categories.append(QgsRendererCategory(6, QgsLineSymbol.createSimple({'color': '#111111', 'width': '3.0'}), "Center Dot"))
        self._strikes_layer.setRenderer(QgsCategorizedSymbolRenderer('id', categories))

        label_settings = QgsPalLayerSettings()
        label_settings.fieldName = "label"
        text_format = QgsTextFormat()
        text_format.setFont(QFont("Arial"))
        text_format.setSizeUnit(QgsUnitTypes.RenderUnit.RenderPixels) 
        text_format.setSize(10.0)  
        text_format.setColor(QColor("#2e7d32"))
        label_settings.setFormat(text_format)
        label_settings.placement = QgsPalLayerSettings.Placement.Line
        label_settings.placementFlags = QgsPalLayerSettings.LinePlacementFlags.AboveLine
        label_settings.obstacle = False                             
        label_settings.dist = 4.0 
        self._strikes_layer.setLabeling(QgsVectorLayerSimpleLabeling(label_settings))
        self._strikes_layer.setLabelsEnabled(True)

        self._strikes_layer.updateExtents()
        self._add_canvas_layer(self._strikes_layer)
        self._update_canvas_visibility()
        
        if hasattr(self, 'iface') and self.iface:
            self.iface.messageBar().pushMessage("Success", f"Projected {len(self.projected_structural_dips)} structural points.", level=0, duration=4)

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
        except Exception: # nosec B110
            pass  # Silently ignore: color extraction failed, return default
        return default_color

    def _compute_projected_traces_nodes(self):
        if not self.dem_layer or not self.trace_layer or not self._profile_layer:
            QMessageBox.warning(self, "No Profile", "Please extract the topographic profile first to initialize data memory.")
            return

        pt4_layer_id = self.combo_pt4_traces_layer.currentData()
        pt4_layer = QgsProject.instance().mapLayer(pt4_layer_id) if pt4_layer_id else None
        if not pt4_layer:
            QMessageBox.warning(self, "Missing Layer", "Please select a valid vector Traces layer.")
            return

        if self.chk_pt4_selected_only.isChecked():
            features = list(pt4_layer.selectedFeatures())
            if not features:
                QMessageBox.warning(self, "No Selection", "No features selected in the chosen layer.")
                return
        else:
            features = list(pt4_layer.getFeatures())

        section_features = list(self.trace_layer.getFeatures())
        if not section_features: return

        section_geom = QgsGeometry(section_features[0].geometry())
        section_geom.get().dropZValue()
        section_geom.get().dropMValue()

        xform_section_to_dem = QgsCoordinateTransform(self.trace_layer.crs(), self.dem_layer.crs(), QgsProject.instance())
        section_geom_dem_crs = QgsGeometry(section_geom)
        section_geom_dem_crs.transform(xform_section_to_dem)

        if section_geom_dem_crs.isMultipart():
            parts = section_geom_dem_crs.asMultiPolyline()
            if parts:
                polyline = parts[0]
            else:
                polyline = []
        else:
            polyline = section_geom_dem_crs.asPolyline()
        if len(polyline) < 2: return

        buffer_dist = self.spin_pt4_buffer.value()
        proj_trend = self.spin_pt4_trend.value()
        proj_plunge = self.spin_pt4_plunge.value()

        xform_traces_to_dem = QgsCoordinateTransform(pt4_layer.crs(), self.dem_layer.crs(), QgsProject.instance())
        dp = self.dem_layer.dataProvider()
        section_length = section_geom_dem_crs.length()

        p_start = polyline[0]
        p_end = polyline[-1]
        A = p_end.y() - p_start.y()
        B = p_start.x() - p_end.x()
        C = p_end.x() * p_start.y() - p_start.x() * p_end.y()

        alpha_proj = math.radians(proj_trend)
        beta_proj = math.radians(proj_plunge)
        px = math.sin(alpha_proj) * math.cos(beta_proj)
        py = math.cos(alpha_proj) * math.cos(beta_proj)
        pz = -math.sin(beta_proj)

        denom = A * px + B * py
        if math.isclose(denom, 0.0):
            QMessageBox.critical(self, "Math Error", "The projection vector is parallel to the section line.")
            return

        chosen_label_field = self.combo_pt4_label_field.currentData()
        projected_lines_features = []

        crs_canvas = self.map_canvas.mapSettings().destinationCrs()
        temp_layer = QgsVectorLayer(f"LineString?crs={crs_canvas.authid()}", "Temp", "memory")
        temp_prov = temp_layer.dataProvider()
        temp_prov.addAttributes([
            QgsField("id", QMetaType.Type.Int), 
            QgsField("label", QMetaType.Type.QString),
            QgsField("color", QMetaType.Type.QString)
        ])
        temp_layer.updateFields()

        for feature in features:
            feat_color_hex = self._get_feature_color_safely(pt4_layer, feature)

            if chosen_label_field and feature[chosen_label_field] is not None:
                lbl_val = str(feature[chosen_label_field])
            else:
                lbl_val = ""

            t_geom = feature.geometry()
            if t_geom.isEmpty(): continue

            for part in t_geom.parts():
                vertices_list = list(part.vertices())
                if len(vertices_list) < 2: continue

                projected_points_in_part = []
                for i in range(len(vertices_list) - 1):
                    p1_raw = QgsPointXY(vertices_list[i])
                    p2_raw = QgsPointXY(vertices_list[i+1])

                    pt1_geom = QgsGeometry.fromPointXY(p1_raw)
                    pt1_geom.transform(xform_traces_to_dem)
                    p1 = pt1_geom.asPoint()

                    pt2_geom = QgsGeometry.fromPointXY(p2_raw)
                    pt2_geom.transform(xform_traces_to_dem)
                    p2 = pt2_geom.asPoint()

                    dist_segment = math.sqrt((p2.x() - p1.x())**2 + (p2.y() - p1.y())**2)
                    if dist_segment < 1e-6: continue

                    steps = 5
                    for step in range(steps + 1):
                        t_ratio = step / float(steps)
                        x_interp = p1.x() + t_ratio * (p2.x() - p1.x())
                        y_interp = p1.y() + t_ratio * (p2.y() - p1.y())
                        pt_curr = QgsPointXY(x_interp, y_interp)

                        val_init, ok_init = dp.sample(pt_curr, 1)
                        z_init = val_init if (ok_init and not math.isnan(val_init)) else 0.0

                        t_proj = -(A * pt_curr.x() + B * pt_curr.y() + C) / denom
                        if abs(t_proj) > buffer_dist: continue

                        x_proj = pt_curr.x() + t_proj * px
                        y_proj = pt_curr.y() + t_proj * py
                        z_proj = z_init + t_proj * pz

                        proj_pt = QgsPointXY(x_proj, y_proj)
                        dist_local = section_geom_dem_crs.lineLocatePoint(QgsGeometry.fromPointXY(proj_pt))

                        if 0.0 <= dist_local <= section_length:
                            projected_points_in_part.append(QgsPointXY(dist_local, z_proj))

                if len(projected_points_in_part) >= 2:
                    feat_line = QgsFeature(temp_layer.fields())
                    feat_line.setGeometry(QgsGeometry.fromPolylineXY(projected_points_in_part))
                    feat_line.setAttributes([7, lbl_val, feat_color_hex])
                    projected_lines_features.append(feat_line)

        if not projected_lines_features:
            QMessageBox.information(self, "No results", "No segments matched criteria within limits.")
            return

        if self._traces_layer:
            self._remove_canvas_layer(self._traces_layer)
            self._traces_layer = None

        self._traces_layer = QgsVectorLayer(f"LineString?crs={crs_canvas.authid()}", "Projected Traces", "memory")
        provider = self._traces_layer.dataProvider()
        provider.addAttributes([
            QgsField("id", QMetaType.Type.Int), 
            QgsField("label", QMetaType.Type.QString),
            QgsField("color", QMetaType.Type.QString)
        ])
        self._traces_layer.updateFields()
        provider.addFeatures(projected_lines_features)
        
        categories = []
        used_colors = {}
        for f in self._traces_layer.getFeatures():
            c_hex = str(f["color"]) if f["color"] else "#37474f"
            used_colors[c_hex] = True

        for c_hex in used_colors.keys():
            sym = QgsLineSymbol.createSimple({'color': c_hex, 'width': '0.8'})
            categories.append(QgsRendererCategory(c_hex, sym, f"Projected Trace ({c_hex})"))

        self._traces_layer.setRenderer(QgsCategorizedSymbolRenderer('color', categories))

        label_settings = QgsPalLayerSettings()
        label_settings.fieldName = "label"
        text_format = QgsTextFormat()
        text_format.setFont(QFont("Arial"))
        text_format.setSizeUnit(QgsUnitTypes.RenderUnit.RenderPixels) 
        text_format.setSize(10.0)  
        text_format.setColor(QColor("#333333"))
        label_settings.setFormat(text_format)
        label_settings.placement = QgsPalLayerSettings.Placement.Line
        label_settings.placementFlags = QgsPalLayerSettings.LinePlacementFlags.AboveLine
        label_settings.obstacle = False                             
        label_settings.dist = 2.0 
        self._traces_layer.setLabeling(QgsVectorLayerSimpleLabeling(label_settings))
        self._traces_layer.setLabelsEnabled(True)

        self._traces_layer.updateExtents()
        self._add_canvas_layer(self._traces_layer)
        self._update_canvas_visibility()
        
        QMessageBox.information(self, "Success", f"Projected {len(projected_lines_features)} segment(s) with custom labels.")

    def _update_profile_layer_renderer(self):
        if not self._profile_layer: return
        categories = []
        
        categories.append(QgsRendererCategory("1", QgsLineSymbol.createSimple({'color': '#000000', 'width': '1.8'}), "Topography"))
        categories.append(QgsRendererCategory("2", QgsLineSymbol.createSimple({'color': '#333333', 'width': '1.0'}), "Axes"))
        categories.append(QgsRendererCategory("3", QgsLineSymbol.createSimple({'color': '#111111', 'width': '1.2'}), "Ticks"))

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
                self.map_canvas.unsetMapTool(self._select_proj_tool)
                self._select_proj_tool = None
            if self._traces_layer:
                self._traces_layer.removeSelection()
                self.map_canvas.refresh()
            self.btn_pt4_select_single.setText("Select & Delete Single Projected Line")

    def _delete_selected_projected_line(self):
        if not self._traces_layer: return
        selected_ids = self._traces_layer.selectedFeatureIds()
        if not selected_ids:
            QMessageBox.information(self, "No Selection", "Please click on a projected line first to select it.")
            return

        to_delete = []
        for fid in selected_ids:
            feat = self._traces_layer.getFeature(fid)
            if feat["id"] == 7:
                to_delete.append(fid)

        if to_delete:
            self._traces_layer.dataProvider().deleteFeatures(to_delete)
            self._traces_layer.removeSelection()
            self._traces_layer.updateExtents()
            self._traces_layer.triggerRepaint()
            self.map_canvas.refresh()
            QMessageBox.information(self, "Cleared", f"Deleted {len(to_delete)} selected projected trace(s).")
        else:
            QMessageBox.warning(self, "Invalid Selection", "Please select a projected trace line inside the canvas.")

    def _clear_projected_traces_nodes(self):
        if not self._traces_layer: return
        provider = self._traces_layer.dataProvider()
        old_ids = [f.id() for f in self._traces_layer.getFeatures() if f["id"] == 7]
        if old_ids:
            provider.deleteFeatures(old_ids)
            self._traces_layer.updateExtents()
            self._traces_layer.triggerRepaint()
            self.map_canvas.refresh()
            QMessageBox.information(self, "Cleared", f"Cleared {len(old_ids)} projected trace segment(s).")
        else:
            QMessageBox.information(self, "Cleared", "No projected traces found to delete.")

    def _setup_tab_folding(self):
        side = QVBoxLayout(self.tab_folding_digitization)

        self.chk_single_element = QCheckBox("Single layer only (Fault/Topography)")
        self.chk_single_element.toggled.connect(self._on_single_element_toggled)
        side.addWidget(self.chk_single_element)
        
        side.addWidget(self.linea_separatrice_folding)
        
        side.addWidget(QLabel("Layer thicknesses:"))
        self.table = QTableWidget(0, 2)
        self.table.setHorizontalHeaderLabels(["Layer", "Thickness"])
        self.table.setColumnWidth(0, 110)
        self.table.setColumnWidth(1, 120)
        side.addWidget(self.table)
        self.table.cellChanged.connect(self._on_thickness_changed)

        self.btn_row = QHBoxLayout()
        self.btn_add = QPushButton("+ layer")
        self.btn_remove = QPushButton("- layer")
        self.btn_row.addWidget(self.btn_add)
        self.btn_row.addWidget(self.btn_remove)
        side.addLayout(self.btn_row)
        self.btn_add.clicked.connect(self._add_layer)
        self.btn_remove.clicked.connect(self._remove_layer)

        side.addWidget(QLabel("Boundary being digitized:"))
        self.combo_digitized = QComboBox()
        side.addWidget(self.combo_digitized)
        self.combo_digitized.currentIndexChanged.connect(self._on_digitized_changed)

        self.btn_start = QPushButton("Start digitizing")
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
        """Metodo dedicato al setup degli strumenti di salvataggio, caricamento ed esportazione."""
        layout = QVBoxLayout(self.tab_import_export)

        # GRUPPO GEOPACKAGE
        self.gpkg_group = QGroupBox("GeoPackage Storage")
        gpkg_layout = QVBoxLayout()
        
        self.btn_export_gpkg = QPushButton("Save All Layers to GeoPackage...")
        self.btn_export_gpkg.setStyleSheet("background-color: #00796b; color: white; font-weight: bold;")
        self.btn_export_gpkg.clicked.connect(self._export_all_to_geopackage)
        gpkg_layout.addWidget(self.btn_export_gpkg)

        self.btn_import_gpkg = QPushButton("Load Layers from GeoPackage...")
        self.btn_import_gpkg.setStyleSheet("background-color: #004d40; color: white; font-weight: bold;")
        self.btn_import_gpkg.clicked.connect(self._import_from_geopackage)
        gpkg_layout.addWidget(self.btn_import_gpkg)
        
        self.gpkg_group.setLayout(gpkg_layout)
        layout.addWidget(self.gpkg_group)

        # LINEA DI SEPARAZIONE
        linea_sep = QFrame()
        linea_sep.setFrameShape(QFrame.Shape.HLine)
        linea_sep.setFrameShadow(QFrame.Shadow.Sunken)
        layout.addWidget(linea_sep)

        # SEZIONE ESPORTAZIONE GRAFICA
        lbl_graphics = QLabel("Graphical Export")
        lbl_graphics.setStyleSheet("font-weight: bold;")
        layout.addWidget(lbl_graphics)

        self.btn_export_svg = QPushButton("Export to A4 SVG with Grid...")
        self.btn_export_svg.setStyleSheet("background-color: #37474f; color: white; font-weight: bold;")
        self.btn_export_svg.clicked.connect(self._export_svg)
        layout.addWidget(self.btn_export_svg)

        layout.addStretch()

    def _regenerate_ghost_layer(self):
        if self._ghost_layer:
            self._remove_canvas_layer(self._ghost_layer)
            self._ghost_layer = None
        if not self.topo_points: return

        crs = self.map_canvas.mapSettings().destinationCrs().authid()
        self._ghost_layer = QgsVectorLayer(f"LineString?crs={crs}", "Profile Base Elements (Ghost)", "memory")
        provider = self._ghost_layer.dataProvider()
        
        provider.addAttributes([QgsField("element_type", QMetaType.Type.Int)])
        self._ghost_layer.updateFields()

        feat_topo = QgsFeature(self._ghost_layer.fields())
        feat_topo.setGeometry(QgsGeometry.fromPolylineXY(self.topo_points))
        feat_topo.setAttributes([1])
        provider.addFeatures([feat_topo])
        self._ghost_layer.updateExtents()

        categories = [QgsRendererCategory(1, QgsLineSymbol.createSimple({'color': '#111111', 'width': '1.5'}), "Topography")]
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
        self.table.setEnabled(not checked)
        self.btn_add.setEnabled(not checked)
        self.btn_remove.setEnabled(not checked)
        self.combo_digitized.setEnabled(not checked)
        self._live_update()

    def _add_canvas_layer(self, layer, zoom_to_it=False):
        QgsProject.instance().addMapLayer(layer, addToLegend=False)
        
        if layer not in self._canvas_layers:
            self._canvas_layers.insert(0, layer)
            
        self.map_canvas.setLayers(self._canvas_layers)
        if zoom_to_it: 
            self.map_canvas.setExtent(layer.extent())
        self.map_canvas.refresh()

    def _remove_canvas_layer(self, layer):
        if layer in self._canvas_layers: 
            self._canvas_layers.remove(layer)
        self.map_canvas.setLayers(self._canvas_layers)
        try: 
            QgsProject.instance().removeMapLayer(layer.id())
        except Exception: # nosec B110
            pass  # Best-effort layer removal; proceed even if removal fails
        self.map_canvas.refresh()

    def _refresh_table(self):
        self.table.blockSignals(True)
        self.table.setRowCount(len(self.layers_thickness))
        for i, t in enumerate(self.layers_thickness):
            item = QTableWidgetItem(f"Layer {i+1}")
            item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.table.setItem(i, 0, item)
            self.table.setItem(i, 1, QTableWidgetItem(str(t)))
        self.table.blockSignals(False)

    def _refresh_combo(self):
        self.combo_digitized.blockSignals(True)
        self.combo_digitized.clear()
        for i in range(len(self.layers_thickness) + 1):
            self.combo_digitized.addItem(f"Top Layer {i+1}" if i < len(self.layers_thickness) else f"Bottom Layer {i}")
        self.combo_digitized.setCurrentIndex(self.digitized_index)
        self.combo_digitized.blockSignals(False)

    def _add_layer(self):
        self.layers_thickness.append(20.0)
        self._refresh_table(); self._refresh_combo(); self._live_update()

    def _remove_layer(self):
        r = self.table.currentRow()
        if r < 0: r = len(self.layers_thickness) - 1
        if len(self.layers_thickness) > 1:
            del self.layers_thickness[r]
            self._refresh_table(); self._refresh_combo(); self._live_update()

    def _on_thickness_changed(self, row, col):
        if col != 1: return
        try:
            v = float(self.table.item(row, col).text().replace(",", "."))
            if v > 0: self.layers_thickness[row] = v
        except (ValueError, AttributeError): pass
        self._live_update()

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
                self.btn_pt4_select_single.setChecked(False)
                self._toggle_select_proj_tool(False)
            self.tool = ParallelFoldTool(self.map_canvas, self)
            self.map_canvas.setMapTool(self.tool)
            self.btn_start.setText("Stop digitalizing")
        else:
            if self.tool: self.tool.reset(); self.map_canvas.unsetMapTool(self.tool); self.tool = None
            self.btn_start.setText("Start digitalizing")

    def _apply_layer_symbology(self):
        if not self._output_layer: return
        categories = []
        for i in range(len(self.layers_thickness) + 1):
            s = QgsLineSymbol.createSimple({'color': self.get_boundary_color(i).name(), 'width': '0.7'})
            categories.append(QgsRendererCategory(int(i), s, f"Layer {i}"))
        s_f = QgsLineSymbol.createSimple({'color': '#141414', 'width': '1.2'})
        categories.append(QgsRendererCategory(999, s_f, "Fault"))
        self._output_layer.setRenderer(QgsCategorizedSymbolRenderer('order', categories))
        self._output_layer.triggerRepaint()

    def _save_layers(self):
        if not self.tool or len(self.tool.vertices) < 2: return
        
        self._block_combo_updates = True
        
        try:
            if not self._output_layer:
                crs = self.map_canvas.mapSettings().destinationCrs().authid()
                self._output_layer = QgsVectorLayer(f"LineString?crs={crs}", "Parallel & Structural Layers", "memory")
                
                self._output_layer.dataProvider().addAttributes([
                    QgsField("distance", QMetaType.Type.Double), 
                    QgsField("order", QMetaType.Type.Int), 
                    QgsField("color", QMetaType.Type.QString)
                ])
                self._output_layer.updateFields()
                self._add_canvas_layer(self._output_layer)

            feats = []
            for i, d, geom in self.tool.last_boundaries:
                f = QgsFeature(self._output_layer.fields())
                f.setGeometry(geom)
                f.setAttributes([0.0, 999 if self.chk_single_element.isChecked() else i, self.get_boundary_color(i).name()])
                feats.append(f)
            self._output_layer.dataProvider().addFeatures(feats)
            self._output_layer.updateExtents()
            self._apply_layer_symbology()
            
            self.map_canvas.refresh()
            self.btn_start.setChecked(False)

        finally:
            self._block_combo_updates = False
            self._populate_layers_combos()

    # =====================================================================
    # SALVATAGGIO / CARICAMENTO DA GEOPACKAGE
    # =====================================================================
    def _export_all_to_geopackage(self):
        layers_to_export = {}
        if self._profile_layer and self._profile_layer.featureCount() > 0:
            layers_to_export["Profile_Layer"] = self._profile_layer
        if self._strikes_layer and self._strikes_layer.featureCount() > 0:
            layers_to_export["Strikes_Layer"] = self._strikes_layer
        if self._intersections_layer and self._intersections_layer.featureCount() > 0:
            layers_to_export["Intersections_Layer"] = self._intersections_layer
        if self._traces_layer and self._traces_layer.featureCount() > 0:
            layers_to_export["Traces_Layer"] = self._traces_layer
        if self._output_layer and self._output_layer.featureCount() > 0:
            layers_to_export["Parallel_Folds_Layer"] = self._output_layer

        if not layers_to_export:
            QMessageBox.warning(self, "No Layers", "There are no layers currently active in the section workspace to save.")
            return

        path, _ = QFileDialog.getSaveFileName(self, "Save Workspace as GeoPackage", "section_workspace.gpkg", "GeoPackage (*.gpkg)")
        if not path:
            return

        if os.path.exists(path):
            try:
                os.remove(path)
            except Exception as e:
                QMessageBox.warning(self, "File Access Error", f"Could not clear existing file: {e}")

        errors = []
        is_first_layer = True
        
        transform_context = QgsProject.instance().transformContext()

        for name, layer in layers_to_export.items():
            options = QgsVectorFileWriter.SaveVectorOptions()
            options.fileEncoding = "UTF-8"
            options.layerName = name
            options.datasourceOptions = ["FORCE_FID=YES"]

            if is_first_layer:
                options.actionOnExistingFile = QgsVectorFileWriter.ActionOnExistingFile.CreateOrOverwriteFile
                is_first_layer = False
            else:
                options.actionOnExistingFile = QgsVectorFileWriter.ActionOnExistingFile.CreateOrOverwriteLayer

            err, err_msg, _, _ = QgsVectorFileWriter.writeAsVectorFormatV3(
                layer,
                path,
                transform_context,
                options
            )

            if err == QgsVectorFileWriter.WriterError.NoError:
                try:
                    style_msg = ""
                    from qgis.PyQt.QtXml import QDomDocument
                    doc = QDomDocument()
                    
                    layer.exportNamedStyle(doc, QgsMapLayer.AllStyleCategories)
                    style_qml = doc.toString()
                    
                    if style_qml:
                        layer.saveStyleToDatabase(
                            name,
                            "Default",
                            "Workspace complete style with advanced labeling",
                            style_qml,
                            True,
                            style_msg
                        )
                except Exception as e:
                    print(f"Non è stato possibile salvare lo stile completo per {name}: {e}")
            else:
                errors.append(f"Layer '{name}' failed: {err_msg} (Error Code: {err})")

        if errors:
            QMessageBox.critical(self, "Export Issues", "Some layers could not be saved:\n" + "\n".join(errors))
        else:
            QMessageBox.information(self, "Success", f"All workspace layers ({len(layers_to_export)}) and complete styles saved successfully to:\n{path}")
            
    def _import_from_geopackage(self):
        path, _ = QFileDialog.getOpenFileName(self, "Load Workspace from GeoPackage", "", "GeoPackage (*.gpkg)")
        if not path:
            return

        gpkg_layers = {}
        import sqlite3
        try:
            conn = sqlite3.connect(path)
            cursor = conn.cursor()
            cursor.execute("SELECT table_name FROM gpkg_contents WHERE data_type='features'")
            rows = cursor.fetchall()
            for r in rows:
                layer_name = r[0]
                gpkg_layers[layer_name] = f"{path}|layername={layer_name}"
            conn.close()
        except Exception as e:
            QMessageBox.critical(self, "Read Error", f"Unable to read file metadata:\n{e}")
            return

        loaded_count = 0
        self._block_combo_updates = True

        try:
            # 1. Ripristino del Profilo
            if "Profile_Layer" in gpkg_layers:
                if self._profile_layer:
                    self._remove_canvas_layer(self._profile_layer)
                self._profile_layer = QgsVectorLayer(gpkg_layers["Profile_Layer"], "Topographic Profile with Grid", "ogr")
                if self._profile_layer.isValid():
                    self.topo_points = []
                    for f in self._profile_layer.getFeatures():
                        if f["id"] == 1:
                            self.topo_points = [QgsPointXY(pt) for pt in f.geometry().vertices()]
                            break
                    if self.topo_points:
                        self.profile_length = self.topo_points[-1].x() - self.topo_points[0].x()
                        self.profile_z_min = min(pt.y() for pt in self.topo_points)
                        self.profile_z_max = max(pt.y() for pt in self.topo_points)

                    self._update_profile_layer_renderer()
                    
                    label_settings = QgsPalLayerSettings()
                    label_settings.fieldName = "label"
                    text_format = QgsTextFormat()
                    text_format.setFont(QFont("Arial"))
                    text_format.setSizeUnit(QgsUnitTypes.RenderUnit.RenderPixels) 
                    text_format.setSize(12.0)  
                    text_format.setColor(QColor("#333333"))
                    label_settings.setFormat(text_format)
                    label_settings.placement = QgsPalLayerSettings.Placement.Line
                    label_settings.placementFlags = QgsPalLayerSettings.LinePlacementFlags.BelowLine
                    label_settings.obstacle = False                             
                    label_settings.allowDegradedPlacement = True                
                    label_settings.priority = 10                                
                    label_settings.dist = 2.0 
                    self._profile_layer.setLabeling(QgsVectorLayerSimpleLabeling(label_settings))
                    self._profile_layer.setLabelsEnabled(True)
                    
                    self._add_canvas_layer(self._profile_layer, zoom_to_it=True)
                    self._regenerate_ghost_layer()
                    loaded_count += 1

            # 2. Ripristino dei Bedding
            if "Strikes_Layer" in gpkg_layers:
                if self._strikes_layer:
                    self._remove_canvas_layer(self._strikes_layer)
                self._strikes_layer = QgsVectorLayer(gpkg_layers["Strikes_Layer"], "Projected Bedding", "ogr")
                if self._strikes_layer.isValid():
                    categories = []
                    categories.append(QgsRendererCategory(5, QgsLineSymbol.createSimple({'color': '#2e7d32', 'width': '1.5'}), "Apparent Dip"))
                    categories.append(QgsRendererCategory(6, QgsLineSymbol.createSimple({'color': '#111111', 'width': '3.0'}), "Center Dot"))
                    self._strikes_layer.setRenderer(QgsCategorizedSymbolRenderer('id', categories))
                    
                    label_settings = QgsPalLayerSettings()
                    label_settings.fieldName = "label"
                    text_format = QgsTextFormat()
                    text_format.setFont(QFont("Arial"))
                    text_format.setSizeUnit(QgsUnitTypes.RenderUnit.RenderPixels) 
                    text_format.setSize(10.0)  
                    text_format.setColor(QColor("#2e7d32"))
                    label_settings.setFormat(text_format)
                    label_settings.placement = QgsPalLayerSettings.Placement.Line
                    label_settings.placementFlags = QgsPalLayerSettings.LinePlacementFlags.AboveLine
                    label_settings.obstacle = False                             
                    label_settings.dist = 4.0 
                    self._strikes_layer.setLabeling(QgsVectorLayerSimpleLabeling(label_settings))
                    self._strikes_layer.setLabelsEnabled(True)
                    
                    self._add_canvas_layer(self._strikes_layer)
                    loaded_count += 1

            # 3. Ripristino delle Intersezioni
            if "Intersections_Layer" in gpkg_layers:
                if self._intersections_layer:
                    self._remove_canvas_layer(self._intersections_layer)
                self._intersections_layer = QgsVectorLayer(gpkg_layers["Intersections_Layer"], "Traces Intersection", "ogr")
                if self._intersections_layer.isValid():
                    sym = QgsLineSymbol.createSimple({'color': '#d32f2f', 'width': '1.0', 'line_style': 'dash'})
                    self._intersections_layer.setRenderer(QgsSingleSymbolRenderer(sym))
                    
                    label_settings = QgsPalLayerSettings()
                    label_settings.fieldName = "label"
                    text_format = QgsTextFormat()
                    text_format.setFont(QFont("Arial"))
                    text_format.setSizeUnit(QgsUnitTypes.RenderUnit.RenderPixels) 
                    text_format.setSize(12.0)  
                    text_format.setColor(QColor("#d32f2f"))
                    label_settings.setFormat(text_format)
                    label_settings.placement = QgsPalLayerSettings.Placement.Line
                    label_settings.placementFlags = QgsPalLayerSettings.LinePlacementFlags.AboveLine
                    label_settings.obstacle = False                             
                    label_settings.allowDegradedPlacement = True                
                    label_settings.priority = 8                                
                    label_settings.dist = 2.0 
                    self._intersections_layer.setLabeling(QgsVectorLayerSimpleLabeling(label_settings))
                    self._intersections_layer.setLabelsEnabled(True)
                    
                    self._add_canvas_layer(self._intersections_layer)
                    loaded_count += 1

            # 4. Ripristino delle Tracce
            if "Traces_Layer" in gpkg_layers:
                if self._traces_layer:
                    self._remove_canvas_layer(self._traces_layer)
                self._traces_layer = QgsVectorLayer(gpkg_layers["Traces_Layer"], "Projected Traces", "ogr")
                if self._traces_layer.isValid():
                    categories = []
                    used_colors = {}
                    for f in self._traces_layer.getFeatures():
                        c_hex = str(f["color"]) if f["color"] else "#37474f"
                        used_colors[c_hex] = True
                    for c_hex in used_colors.keys():
                        sym = QgsLineSymbol.createSimple({'color': c_hex, 'width': '0.8'})
                        categories.append(QgsRendererCategory(c_hex, sym, f"Projected Trace ({c_hex})"))
                    self._traces_layer.setRenderer(QgsCategorizedSymbolRenderer('color', categories))
                    
                    label_settings = QgsPalLayerSettings()
                    label_settings.fieldName = "label"
                    text_format = QgsTextFormat()
                    text_format.setFont(QFont("Arial"))
                    text_format.setSizeUnit(QgsUnitTypes.RenderUnit.RenderPixels) 
                    text_format.setSize(10.0)  
                    text_format.setColor(QColor("#333333"))
                    label_settings.setFormat(text_format)
                    label_settings.placement = QgsPalLayerSettings.Placement.Line
                    label_settings.placementFlags = QgsPalLayerSettings.LinePlacementFlags.AboveLine
                    label_settings.obstacle = False                             
                    label_settings.dist = 2.0 
                    self._traces_layer.setLabeling(QgsVectorLayerSimpleLabeling(label_settings))
                    self._traces_layer.setLabelsEnabled(True)
                    
                    self._add_canvas_layer(self._traces_layer)
                    loaded_count += 1

            # 5. Ripristino del Layer Pieghe Parallele (Disegnato)
            if "Parallel_Folds_Layer" in gpkg_layers:
                if self._output_layer:
                    self._remove_canvas_layer(self._output_layer)
                
                temp_ogr_layer = QgsVectorLayer(gpkg_layers["Parallel_Folds_Layer"], "temp_ogr", "ogr")
                if temp_ogr_layer.isValid():
                    crs = temp_ogr_layer.crs().authid()
                    self._output_layer = QgsVectorLayer(f"LineString?crs={crs}", "Parallel & Structural Layers", "memory")
                    provider = self._output_layer.dataProvider()
                    
                    from qgis.core import QgsField, QgsFeature
                    
                    provider.addAttributes([
                        QgsField("type", QMetaType.Type.QString),
                        QgsField("order", QMetaType.Type.Int),
                        QgsField("color", QMetaType.Type.QString)
                    ])
                    self._output_layer.updateFields()
                    
                    fields = self._output_layer.fields()
                    old_fields = temp_ogr_layer.fields()
                    idx_type = old_fields.indexFromName("type")
                    idx_order = old_fields.indexFromName("order")
                    idx_color = old_fields.indexFromName("color")
                    
                    new_features = []
                    for old_feat in temp_ogr_layer.getFeatures():
                        new_feat = QgsFeature(fields)
                        new_feat.setGeometry(old_feat.geometry())
                        
                        new_feat["type"] = old_feat.attribute(idx_type) if idx_type != -1 else ""
                        new_feat["order"] = old_feat.attribute(idx_order) if idx_order != -1 else 0
                        new_feat["color"] = old_feat.attribute(idx_color) if idx_color != -1 else ""
                        
                        new_features.append(new_feat)
                    
                    provider.addFeatures(new_features)
                    self._output_layer.updateExtents()
                    
                    self._apply_layer_symbology()
                    self._add_canvas_layer(self._output_layer)
                    loaded_count += 1
                    

            self._update_canvas_visibility()
            self.map_canvas.refresh()

        finally:
            self._block_combo_updates = False
            self._populate_layers_combos()

        if loaded_count > 0:
            QMessageBox.information(self, "Success", f"Successfully loaded {loaded_count} layers from the GeoPackage!")
        else:
            QMessageBox.warning(self, "No Valid Layers", "No workspace-compatible layers found inside the selected GeoPackage.")

    # =====================================================================

    def _toggle_edit_tool(self, checked):
        if checked:
            if self.btn_start.isChecked(): self.btn_start.setChecked(False)
            if hasattr(self, 'btn_pt4_select_single') and self.btn_pt4_select_single.isChecked():
                self.btn_pt4_select_single.setChecked(False)
                self._toggle_select_proj_tool(False)
            self._edit_tool = AdvancedEditSavedLineTool(self.map_canvas, self)
            self.map_canvas.setMapTool(self._edit_tool)
        else:
            if self._edit_tool: 
                self._edit_tool.clear_markers()
                self.map_canvas.unsetMapTool(self._edit_tool)
                self._edit_tool = None
            if self._output_layer:
                self._output_layer.removeSelection()
                self._output_layer.triggerRepaint()
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
        if ids:
            self._output_layer.dataProvider().deleteFeatures(ids)
            self._output_layer.updateExtents(); self._output_layer.removeSelection()
            if self._edit_tool: self._edit_tool.clear_markers()
            self.map_canvas.refresh()

    @staticmethod
    def _calculate_nice_step(delta):
        raw = delta / 6.0
        if raw == 0: return 1.0
        log = math.floor(math.log10(raw))
        b = 10 ** log
        rel = raw / b
        nice = 1.0 if rel < 1.5 else 2.0 if rel < 3.5 else 5.0 if rel < 7.5 else 10.0
        return nice * b

    def _export_svg(self):
        data_to_export = []
        texts_to_export = []  
        xs, ys = [], []

        if self._output_layer and self._output_layer.featureCount() > 0:
            for f in self._output_layer.getFeatures():
                if f.geometry().isEmpty(): continue
                pts = [QPointF(p.x(), p.y()) for p in f.geometry().vertices()]
                for p in f.geometry().vertices():
                    xs.append(p.x()); ys.append(p.y())
                data_to_export.append((pts, QColor(f["color"]), 0.5))

        if self._profile_layer and self._profile_layer.featureCount() > 0:
            for f in self._profile_layer.getFeatures():
                if f.geometry().isEmpty(): continue
                geom = f.geometry()
                pts = [QPointF(p.x(), p.y()) for p in geom.vertices()]
                for p in geom.vertices():
                    xs.append(p.x()); ys.append(p.y())
                
                prof_id = f["id"]
                color = QColor("#000000")
                width = 0.4
                
                if prof_id == 1: 
                    color = QColor("#000000")
                    width = 0.7
                elif prof_id == 2: 
                    color = QColor("#333333")
                    width = 0.3
                elif prof_id == 3: 
                    color = QColor("#555555")
                    width = 0.2
                    label = str(f["label"]) if f["label"] is not None else ""
                    if label:
                        vertices = list(geom.vertices())
                        if len(vertices) >= 2:
                            p_text = vertices[1]
                            is_y_axis = (vertices[0].x() == vertices[1].x())
                            align = "right" if is_y_axis else "center"
                            offset_x = -1.5 if is_y_axis else 0.0
                            offset_y = 0.0 if is_y_axis else -2.5
                            texts_to_export.append((
                                label, 
                                p_text.x() + offset_x, 
                                p_text.y() + offset_y, 
                                QColor("#333333"), 
                                3.0, 
                                align
                            ))

                data_to_export.append((pts, color, width))

        if self._intersections_layer and self._intersections_layer.featureCount() > 0:
            for f in self._intersections_layer.getFeatures():
                if f.geometry().isEmpty(): continue
                geom = f.geometry()
                pts = [QPointF(p.x(), p.y()) for p in geom.vertices()]
                for p in geom.vertices():
                    xs.append(p.x()); ys.append(p.y())
                data_to_export.append((pts, QColor("#d32f2f"), 0.4))

                label = str(f["label"]) if f["label"] is not None else ""
                if label:
                    vertices = list(geom.vertices())
                    if vertices:
                        top_pt = max(vertices, key=lambda p: p.y())
                        texts_to_export.append((
                            label, 
                            top_pt.x(), 
                            top_pt.y() + 2.0, 
                            QColor("#d32f2f"), 
                            3.5, 
                            "center"
                        ))

        if self._strikes_layer and self._strikes_layer.featureCount() > 0:
            for f in self._strikes_layer.getFeatures():
                if f.geometry().isEmpty(): continue
                geom = f.geometry()
                pts = [QPointF(p.x(), p.y()) for p in geom.vertices()]
                for p in geom.vertices():
                    xs.append(p.x()); ys.append(p.y())
                col = QColor("#2e7d32") if f["id"] == 5 else QColor("#111111")
                width = 0.5 if f["id"] == 5 else 0.8
                data_to_export.append((pts, col, width))

                label = str(f["label"]) if f["label"] is not None else ""
                if label and f["id"] == 5:
                    vertices = list(geom.vertices())
                    if vertices:
                        mid_x = sum(p.x() for p in vertices) / len(vertices)
                        mid_y = sum(p.y() for p in vertices) / len(vertices)
                        texts_to_export.append((
                            label, 
                            mid_x, 
                            mid_y + 3.0, 
                            QColor("#2e7d32"), 
                            2.8, 
                            "center"
                        ))

        if self._traces_layer and self._traces_layer.featureCount() > 0:
            for f in self._traces_layer.getFeatures():
                if f.geometry().isEmpty(): continue
                geom = f.geometry()
                pts = [QPointF(p.x(), p.y()) for p in geom.vertices()]
                for p in geom.vertices():
                    xs.append(p.x()); ys.append(p.y())
                color = QColor(f["color"]) if f["color"] else QColor("#37474f")
                data_to_export.append((pts, color, 0.3))

                label = str(f["label"]) if f["label"] is not None else ""
                if label:
                    vertices = list(geom.vertices())
                    if vertices:
                        mid_pt = vertices[len(vertices) // 2]
                        texts_to_export.append((
                            label, 
                            mid_pt.x(), 
                            mid_pt.y() + 2.0, 
                            QColor("#333333"), 
                            2.8, 
                            "center"
                        ))

        if not xs:
            QMessageBox.warning(self, "Export Error", "No vector elements found to export.")
            return

        path, _ = QFileDialog.getSaveFileName(self, "Export SVG", "section_complete.svg", "SVG (*.svg)")
        if not path: return

        pw, ph, mar = 297.0, 210.0, 25.0
        x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
        
        sc = min((pw - 2*mar)/(x1-x0 or 1), (ph - 2*mar-10)/(y1-y0 or 1))
        ox = mar + ((pw-2*mar)-(x1-x0)*sc)/2
        oy = mar + ((ph-2*mar-10)-(y1-y0)*sc)/2 + 10

        def to_svg(x, y): 
            return QPointF(ox + (x-x0)*sc, (ph-oy) - (y-y0)*sc)

        gen = QSvgGenerator()
        gen.setFileName(path)
        gen.setSize(QSize(int(pw*3.77), int(ph*3.77)))
        gen.setViewBox(QRectF(0, 0, pw, ph))
        
        p = QPainter()
        if p.begin(gen):
            p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
            p.fillRect(QRectF(0, 0, pw, ph), QBrush(Qt.GlobalColor.white))
            
            p.setPen(QPen(QColor(180, 180, 180), 0.15))
            p.drawRect(QRectF(to_svg(x0, y1), to_svg(x1, y0)))
            
            for pts, col, width in data_to_export:
                p.setPen(QPen(col, width))
                p.drawPolyline(QPolygonF([to_svg(pt.x(), pt.y()) for pt in pts]))
            
            for text, x_map, y_map, col, size_mm, align in texts_to_export:
                pt_canvas = to_svg(x_map, y_map)
                
                font = QFont("Arial")
                font.setPointSizeF(size_mm) 
                p.setFont(font)
                p.setPen(QPen(col, 0.1))
                
                fm = p.fontMetrics()
                tw = fm.horizontalAdvance(text)
                th = fm.height()
                
                if align == "center":
                    text_x = pt_canvas.x() - (tw / 2.0)
                elif align == "right":
                    text_x = pt_canvas.x() - tw - 1.0
                else: 
                    text_x = pt_canvas.x() + 1.0
                
                text_y = pt_canvas.y() + (th / 3.0) 
                
                p.drawText(QPointF(text_x, text_y), text)
            
            p.end()
            QMessageBox.information(self, "Success", "Complete profile section exported correctly to SVG.")
            
    def _clear_points(self):
        if self.tool:
            self.tool.clear_points()

    def closeEvent(self, event):
        if self.tool: self.map_canvas.unsetMapTool(self.tool)
        if self._edit_tool: self._edit_tool.clear_markers(); self.map_canvas.unsetMapTool(self._edit_tool)
        if self._select_proj_tool: self.map_canvas.unsetMapTool(self._select_proj_tool)
        for l in list(self._canvas_layers):
            try: QgsProject.instance().removeMapLayer(l.id())
            except: pass
        super().closeEvent(event)


class ParallelFoldTool(QgsMapTool):
    def __init__(self, canvas, window):
        super().__init__(canvas)
        self.canvas, self.window = canvas, window
        self.vertices, self.vertex_markers, self.rubber_bands, self.last_boundaries = [], [], [], []
        self.dragged_index = None

    def clear_points(self):
        self.vertices = []
        for m in self.vertex_markers: 
            self.canvas.scene().removeItem(m)
        self.vertex_markers = []
        for rb in self.rubber_bands: 
            rb.reset(QgsWkbTypes.GeometryType.LineGeometry)
        self.rubber_bands = []

    def set_vertices(self, nv):
        self.clear_points()
        self.vertices = list(nv)
        for pt in self.vertices:
            m = QgsVertexMarker(self.canvas)
            m.setCenter(pt)
            m.setIconType(QgsVertexMarker.ICON_CIRCLE)
            m.setIconSize(8)
            m.setPenWidth(2)
            m.setColor(QColor(55,138,221))
            self.vertex_markers.append(m)
        self.recalc_and_draw()

    def reverse_vertices(self):
        self.vertices.reverse()
        self.vertex_markers.reverse()
        self.recalc_and_draw()

    def reset(self): 
        self.clear_points()

    def _hit_test(self, pos):
        for i, pt in enumerate(self.vertices):
            sc = self.toCanvasCoordinates(pt)
            if math.hypot(sc.x() - pos.x(), sc.y() - pos.y()) <= 12: 
                return i
        return None

    def canvasPressEvent(self, e):
        if e.button() == Qt.MouseButton.RightButton: 
            return
        
        idx = self._hit_test(e.pos())
        if idx is not None: 
            self.dragged_index = idx
            return
        
        pt = self.toMapCoordinates(e.pos())
        pt_xy = QgsPointXY(pt)
        inserted = False
        
        if len(self.vertices) >= 2:
            geom_line = QgsGeometry.fromPolylineXY(self.vertices)
            tolerance = 12 * self.canvas.mapUnitsPerPixel()
            
            if geom_line.distance(QgsGeometry.fromPointXY(pt_xy)) <= tolerance:
                best_idx = -1
                min_dist = float('inf')
                
                for i in range(len(self.vertices) - 1):
                    p1 = self.vertices[i]
                    p2 = self.vertices[i+1]
                    seg_geom = QgsGeometry.fromPolylineXY([p1, p2])
                    dist = seg_geom.distance(QgsGeometry.fromPointXY(pt_xy))
                    if dist < min_dist:
                        min_dist = dist
                        best_idx = i + 1
                
                if best_idx != -1:
                    self.vertices.insert(best_idx, pt_xy)
                    m = QgsVertexMarker(self.canvas)
                    m.setCenter(pt_xy)
                    m.setIconType(QgsVertexMarker.ICON_CIRCLE)
                    m.setIconSize(8)
                    m.setColor(QColor(55,138,221))
                    self.vertex_markers.insert(best_idx, m)
                    inserted = True

        if not inserted:
            self.vertices.append(pt_xy)
            m = QgsVertexMarker(self.canvas)
            m.setCenter(pt_xy)
            m.setIconType(QgsVertexMarker.ICON_CIRCLE)
            m.setIconSize(8)
            m.setColor(QColor(55,138,221))
            self.vertex_markers.append(m)
            
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
        if e.button() == Qt.MouseButton.RightButton:
            return
            
        pos = e.pos()
        idx = self._hit_test(pos)
        
        if idx is not None:
            if len(self.vertices) > 2:
                self.vertices.pop(idx)
                m_to_remove = self.vertex_markers.pop(idx)
                self.canvas.scene().removeItem(m_to_remove)
                self.recalc_and_draw()

    def recalc_and_draw(self):
        for rb in self.rubber_bands: 
            rb.reset(QgsWkbTypes.GeometryType.LineGeometry)
        self.rubber_bands = []
        self.last_boundaries = []
        
        if len(self.vertices) < 2: 
            return
            
        self.window.update_ghost_bounds_from_digitizing(
            min(p.y() for p in self.vertices), 
            max(p.y() for p in self.vertices)
        )
        offsets = self.window.cumulative_offsets()
        lg = QgsGeometry.fromPolylineXY(self.vertices)
        sign = 1.0 if self.vertices[-1].x() >= self.vertices[0].x() else -1.0
        
        for i, d in enumerate(offsets):
            ad = d * sign
            geom = QgsGeometry(lg) if abs(ad) < 1e-9 else lg.offsetCurve(ad, 8, JOIN_STYLE_ROUND, 2.0)
            if not geom or geom.isEmpty(): 
                continue
            self.last_boundaries.append((i, d, geom))
            rb = QgsRubberBand(self.canvas, QgsWkbTypes.GeometryType.LineGeometry)
            rb.setToGeometry(geom, None)
            rb.setColor(self.window.get_boundary_color(i))
            rb.setWidth(4 if i == self.window.digitized_index else 2)
            self.rubber_bands.append(rb)


class AdvancedEditSavedLineTool(QgsMapTool):
    def __init__(self, canvas, window):
        super().__init__(canvas)
        self.canvas, self.window = canvas, window
        self.feat, self.markers, self.dragged_marker = None, [], None
        self.dragged_idx = None

    def clear_markers(self):
        for m in self.markers:
            if m:
                self.canvas.scene().removeItem(m)
        self.markers = []

    def deactivate(self):
        self.clear_markers()
        if self.window._output_layer:
            self.window._output_layer.removeSelection()
            self.window._output_layer.triggerRepaint()
        self.canvas.refresh()
        super().deactivate()

    def canvasPressEvent(self, e):
        if self.window.main_tabs.currentIndex() != 1:
            return
        layer = self.window._output_layer
        if not layer: return
        pt = self.toMapCoordinates(e.pos())
        
        for i, m in enumerate(self.markers):
            if math.hypot(m.center().x() - pt.x(), m.center().y() - pt.y()) < 15 * self.canvas.mapUnitsPerPixel():
                self.dragged_marker = m
                self.dragged_idx = i
                return

        cg = QgsGeometry.fromPointXY(pt)
        bf, bd = None, 12 * self.canvas.mapUnitsPerPixel()
        for f in layer.getFeatures():
            d = f.geometry().distance(cg)
            if d <= bd: bd, bf = d, f
        
        if bf:
            layer.selectByIds([bf.id()])
            self.feat = bf
            self.clear_markers()
            for p in bf.geometry().vertices():
                m = QgsVertexMarker(self.canvas)
                m.setCenter(QgsPointXY(p))
                m.setIconType(QgsVertexMarker.ICON_BOX)
                m.setColor(QColor(230, 20, 150))
                self.markers.append(m)
        else:
            layer.removeSelection()
            self.feat = None
            self.clear_markers()
        
        layer.triggerRepaint()
        self.canvas.refresh()

    def canvasMoveEvent(self, e):
        if self.dragged_marker:
            pt = self.toMapCoordinates(e.pos())
            self.dragged_marker.setCenter(QgsPointXY(pt))

    def canvasReleaseEvent(self, e):
        if self.dragged_marker and self.feat:
            new_pts = [QgsPointXY(m.center()) for m in self.markers]
            self.window._output_layer.dataProvider().changeGeometryValues({
                self.feat.id(): QgsGeometry.fromPolylineXY(new_pts)
            })
            self.window._output_layer.triggerRepaint()
            self.canvas.refresh()
            self.dragged_marker = None
        
    def canvasDoubleClickEvent(self, e):
        if self.feat and self.dragged_idx is not None:
            if len(self.markers) > 2:
                m_to_remove = self.markers.pop(self.dragged_idx)
                self.canvas.scene().removeItem(m_to_remove)
                
                new_pts = [QgsPointXY(m.center()) for m in self.markers]
                self.window._output_layer.dataProvider().changeGeometryValues({
                    self.feat.id(): QgsGeometry.fromPolylineXY(new_pts)
                })
                self.window._output_layer.triggerRepaint()
                self.canvas.refresh()
            
            self.dragged_marker = None
            self.dragged_idx = None


class SelectProjectedLineTool(QgsMapTool):
    def __init__(self, canvas, window):
        super().__init__(canvas)
        self.canvas = canvas
        self.window = window

    def canvasPressEvent(self, e):
        layer = self.window._traces_layer
        if not layer: return

        pt = self.toMapCoordinates(e.pos())
        cg = QgsGeometry.fromPointXY(pt)
        
        search_distance = 10 * self.canvas.mapUnitsPerPixel()
        
        best_feat = None
        min_dist = search_distance

        for f in layer.getFeatures():
            if f["id"] != 7:
                continue
            
            dist = f.geometry().distance(cg)
            if dist <= min_dist:
                min_dist = dist
                best_feat = f

        if best_feat:
            layer.selectByIds([best_feat.id()])
        else:
            layer.removeSelection()
            
        self.canvas.refresh()

    def keyPressEvent(self, e):
        if e.key() in (Qt.Key.Key_Delete, Qt.Key.Key_Backspace):
            self.window._delete_selected_projected_line()
