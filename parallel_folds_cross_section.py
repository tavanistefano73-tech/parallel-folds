# -*- coding: utf-8 -*-
from qgis.PyQt.QtWidgets import QAction
from qgis.PyQt.QtGui import QIcon
from qgis.core import Qgis


class ParallelFoldsPlugin:
    """
    Entry point del plugin. Rileva automaticamente la versione di QGIS
    in esecuzione e carica l'implementazione compatibile:
      - QGIS >= 4.0  -> core_qgis4.py (PyQt6 / API scoped)
      - QGIS  < 4.0  -> core_qgis3.py (PyQt5 / API "flat")
    """

    def __init__(self, iface):
        self.iface = iface
        self.actions = []
        self.menu = "&Parallel Folds Tool"
        self.dlg = None  # riferimento persistente, altrimenti la finestra
                          # verrebbe distrutta dal garbage collector

    def initGui(self):
        icon = QIcon(":/images/themes/default/mIconPluginsMenu.svg")
        self.action = QAction(icon, "Parallel Folds & Structural Tool", self.iface.mainWindow())
        self.action.triggered.connect(self.run)
        self.iface.addToolBarIcon(self.action)
        self.iface.addPluginToMenu(self.menu, self.action)
        self.actions.append(self.action)

    def unload(self):
        for action in self.actions:
            self.iface.removePluginMenu(self.menu, action)
            self.iface.removeToolBarIcon(action)
        self.actions = []

    def _is_qgis4_or_newer(self):
        # Qgis.QGIS_VERSION_INT ha il formato MNNPP (es. 34099 = 3.40.99, 40200 = 4.2.0)
        return Qgis.QGIS_VERSION_INT >= 40000

    def run(self):
        if self._is_qgis4_or_newer():
            from .core_qgis4 import FoldWindow
        else:
            from .core_qgis3 import FoldWindow

        # Se la finestra e' gia' aperta, la porta in primo piano invece di
        # crearne una seconda istanza
        if self.dlg is not None:
            try:
                self.dlg.show()
                self.dlg.raise_()
                self.dlg.activateWindow()
                return
            except RuntimeError:
                # la finestra precedente e' stata chiusa/distrutta
                self.dlg = None

        self.dlg = FoldWindow(self.iface)
        self.dlg.show()
