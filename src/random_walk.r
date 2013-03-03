##############################################################################
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
##############################################################################

# Run as follows: R --slave --args node_file edge_file output_file < random_walk.r

# Parameter r: restart probability
r<-0.75
# Parameter max_n_iter: max number of iterations for convergence 
max_n_iter<-50
# Parameter r: convergence cutoff 
convergence_cutoff<-1e-6

main<-function() {
    #test_with_toy_data()
    #test_with_file_data()
    args = commandArgs(trailingOnly=T)
    node_file = args[1]
    edge_file = args[2]
    out_file = args[3]
    print(c(node_file, edge_file, out_file))
    if(length(args)==3) {
	run_random_walk_with_restart(node_file, edge_file, out_file)
    } else if(length(args)==4) {
	run_random_walk_with_restart(node_file, edge_file, out_file, propagation=TRUE)
    }
}

run_random_walk_with_restart<-function(node_file, edge_file, out_file, propagation=FALSE) {
    A<-get_adjacency_matrix(edge_file)
    p_0<-get_node_scores(node_file)
    p_tx<-random_walk_with_restart(A, p_0, r, max_n_iter, convergence_cutoff, propagation)
    write.table(p_tx$scores, out_file, col.names=F, quote=F)
}

get_adjacency_matrix<-function(edge_file) {
    #A<-as.matrix(read.table(edge_file, header = T, row.names=1, sep = " ", skip = 0))
    #return(A)
    d<-read.table(edge_file, header = F, sep = " ", col.names=c("id1", "score", "id2"), skip = 0)
    ids<-unique(c(as.vector(d$id1),as.vector(d$id2)))
    ids<-sort(ids)
    n<-length(ids)
    A<-matrix(rep(0, n*n), n, n)
    rownames(A)<-ids
    colnames(A)<-ids
    for (i in 1:nrow(d)) {
	A[as.character(d[i,"id1"]), as.character(d[i,"id2"])] = d[i,"score"]
	A[as.character(d[i,"id2"]), as.character(d[i,"id1"])] = d[i,"score"]
    }
    return(A)
}

get_node_scores<-function(node_file) {
    d<-read.table(node_file, header = F, sep = " ", col.names=c("id", "score"), skip = 0)
    ids<-as.vector(d$id)
    ids<-sort(ids)
    scores<-rep(0, length(ids))
    names(scores)<-ids
    for (i in 1:nrow(d)) {
	scores[as.character(d[i,"id"])] = d[i,"score"]
    }
    return(scores)
}

random_walk_with_restart <- function(A, p_0, r, max_n_iter, convergence_cutoff, propagation=FALSE) {

    if(propagation) {
	# Weighted degrees W(v) = sum(w(u,v)) for u,v in edges
	w<-colSums(A)
    	#w<-rep(0, nrow(A))
	#for (i in 1:nrow(A)) {
	#    w[i]<-sum(A[i,])
	#}
	# Degree weighted adjacency matrix = w(uv) / sqrt(W(u)*W(v))
	W<-matrix(rep(1,nrow(A)*ncol(A)), nrow(A), ncol(A))
	for (i in 1:nrow(A)) {
	    for (j in 1:ncol(A)) {
		v<-sqrt(w[i]*w[j])
		W[i,j]<-v
		W[j,i]<-v
	    }
	}
	W<-A/W
    } else {
	# Convert A to column-normalized adjacency matrix 
	W<-scale(A,center=F,scale=colSums(A))
    }

    # Assign equal probabilities to seed nodes
    p_0<-p_0/sum(p_0)
    n<-length(p_0)
    dim(p_0)<-c(n,1)
    #p_t<-rep(0, n)
    #dim(p_t)<-c(n,1)
    p_t<-p_0

    # Iterate till convergance is met or max_n_iter is exceeded
    for ( i in 1:max_n_iter ) {
	# Calculate new proabalities
	p_tx<- (1-r) * W %*% p_t + r * p_0
	# Check convergance
	if ( norm(p_tx-p_t) < convergence_cutoff ) {
	#if ( norm(p_tx-p_t, "F") < convergence_cutoff ) {
		break
	}
	p_t<-p_tx
    }

    return(list(scores=p_tx, iter=i))
}

test_with_file_data<-function() {

    edge_file = "../data/input_runs_for_draft/entrez/edge_scores.sif"
    node_file = "../data/input_runs_for_draft/entrez/chen_autism/node_scores.sif"
    #edge_file = "test_ppi.dat"

    #edge_file = "../data/toy_data/test_interactions_small.sif"
    #node_file = "../data/toy_data/test_proteins_small.sif"

    A<-get_adjacency_matrix(edge_file)
    print(dim(A))
    p_0<-get_node_scores(node_file)
    print(length(p_0))

    p_tx<-random_walk_with_restart(A, p_0, r, max_n_iter, convergence_cutoff)
    print(p_tx)
}

test_with_toy_data<-function() {

    M<-matrix(rep(0,25),ncol=5)
    M[1,]=c(0,1,0,0,0)
    M[2,]=c(1,0,1,1,0)
    M[3,]=c(0,1,0,1,0)
    M[4,]=c(0,1,1,0,1)
    M[5,]=c(0,0,0,1,0)

    p_0<-t(t(c(0,1,1,0,0)))

    p_tx<-random_walk_with_restart(M, p_0, r, max_n_iter, convergence_cutoff)
    print(p_tx)
}

main()

