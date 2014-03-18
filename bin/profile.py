#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Profiling the GFF parsing program and extracting upstream sequences
import profile

profile.run("from gff import main; main.main()","../profile1803.tmp")
