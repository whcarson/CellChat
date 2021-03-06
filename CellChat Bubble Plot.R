## Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(LCMA01_BEVR_CD31_Chat, sources.use = c(1:12), targets.use = c(1:12), remove.isolate = FALSE)


