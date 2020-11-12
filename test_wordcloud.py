import wordcloud

txt = ['tissue:flowers+ecotype:Col-0', 'ecotype:Col-0+genotype/variation:wildtype+rnatype: polysomeassociatedRNA+librarytype:QuantSeq3\'endmRNA-seq']
txt = '--'.join(txt)
w = wordcloud.WordCloud( width=1000, height=700,background_color="white")
w.generate(txt)
w.to_file("map.png")

