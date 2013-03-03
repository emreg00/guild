#############################################################################
#    GUILD (Genes Underlying Inheritance Linked Disorders) implements several 
#    graph based algorithms for scoring relevance of a node in the network in 
#    terms of a phenotype using known associations in the node's neighborhood 
#    for that phenotype. GUILD has been applied to the prioritization of genes 
#    for several human disorders. 2011 - Emre Guney (Unviersitat Pompeu Fabra)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#############################################################################

import sys
import calculate_mean_and_sigma

def main():
    scores_file_list = sys.argv[1:-1]
    output_scores_file = sys.argv[-1]
    score_combined(scores_file_list, output_scores_file)
    return

def score_combined(scores_file_list, output_scores_file):
    """
	Calculates a combined score based on normalized scores of each scoring method
    """
    node_to_scores = {}
    for scores_file in scores_file_list:
	node_to_score_inner = {}
	for line in open(scores_file):
	    node, score = line.strip().split() 
	    node_to_score_inner[node] = float(score)
	mean, sigma = calculate_mean_and_sigma.calc_mean_and_sigma(node_to_score_inner.values())
	for node, score in node_to_score_inner.iteritems():
	    node_to_scores.setdefault(node, []).append((score-mean)/sigma)
    values = []
    for node, scores in node_to_scores.iteritems():
	score = sum(scores) / len(scores)
	values.append((score, node))
    values.sort()
    min_v, max_v = min(values)[0], max(values)[0]
    f = open(output_scores_file, 'w')
    for score, node in values:
	score = (score-min_v) / (max_v-min_v)
	f.write("%s\t%f\n" % (node, score))
    f.close()
    return


if __name__ == "__main__":
    main()


