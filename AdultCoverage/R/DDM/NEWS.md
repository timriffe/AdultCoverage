
# DDM version 1.2.xxx
- [x] several line fitting choices offered for ggb and ggbseg methods (new arg = lm.method)
- [x] two options for intercensal birthdays per annum (new arg = nx.method)
- [x] coverage calculations corrected, more transparent
- [ ] coverage calculations better documented
- [ ] automatic age trim approach documented
- [ ] age varying coverage methods added

# DDM version 1.1.0
- SEG method now has its own automatic age trim selection
- GGBSEG method now has two age trims for the first and second stages, each with automatic solutions
- ggbChooseAges() ages outside maxA and minA no longer plotted

# DDM version 1.0.1
- change to BH closeout. Now consistent with source spreadsheet. HT Haidong Wang for pointing this out.

# DDM version 1.0.0.0

- first release! methods include ggb(), seg(), and ggbseg(), as well as diagnostic functions ggbChooseAges() and segplot().

