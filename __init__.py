# -*- coding: utf-8 -*-
def classFactory(iface):
    from .parallel_folds_plugin import ParallelFoldsPlugin
    return ParallelFoldsPlugin(iface)
