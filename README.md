# AlignmentFreeD2
Python implementation of D2 associated statistics evaluting sequence similarities

According to *Alignment-Free Sequence Comparison (I):Statistics and Power*, D2* and D2S outperform the original D2 statistic. We implemented D2* and D2S in python.


## Usage:

Users should first build the sequence occurrence probability model with a list of sequences.

Then users can compute any two sequence's similarity using the probability model


## ToDo:
- [ ] Design a probablistic similarity evaluation base class
- [ ] Design a D2*, D2S class based on the base class
- [ ] Implement the classes above
- [ ] Design higher-order Markov sequence models 
- [ ] Implement Markov chain model estimation and evaluation with pomegranate
