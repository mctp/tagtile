tagv3.sel.tile <- tagv3.sel[(edge)]
tagv3.sel.tile.gr <- with(tagv3.sel.tile, GRanges(transcript_id, IRanges(beg+1, end), probe_id=probe_id))
tagv3.sel.tile.gr.onGm <- mapFromTranscripts(tagv3.sel.tile.gr, gen26.sel1.tx)
export(tagv3.sel.tile.gr.onGm, "edge.bed")

tagv3.sel.tile <- tagv3.sel[(grid)]
tagv3.sel.tile.gr <- with(tagv3.sel.tile, GRanges(transcript_id, IRanges(beg+1, end), probe_id=probe_id))
tagv3.sel.tile.gr.onGm <- mapFromTranscripts(tagv3.sel.tile.gr, gen26.sel1.tx)
export(tagv3.sel.tile.gr.onGm, "grid.bed")

tagv3.sel.tile <- tagv3.sel[(cosmic)]
tagv3.sel.tile.gr <- with(tagv3.sel.tile, GRanges(transcript_id, IRanges(beg+1, end), probe_id=probe_id))
tagv3.sel.tile.gr.onGm <- mapFromTranscripts(tagv3.sel.tile.gr, gen26.sel1.tx)
export(tagv3.sel.tile.gr.onGm, "cosmic.bed")

