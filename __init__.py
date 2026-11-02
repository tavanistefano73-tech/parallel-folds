# -*- coding: utf-8 -*-
def classFactory(iface):
    """Entry point required by QGIS."""
    from .parallel_folds_cross_section import ParallelFoldsPlugin
    return ParallelFoldsPlugin(iface)