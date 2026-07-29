# -*- coding: utf-8 -*-
def classFactory(iface):
    from .parallel_folds_cross_section import ParallelFoldsPlugin
    return ParallelFoldsPlugin(iface)
