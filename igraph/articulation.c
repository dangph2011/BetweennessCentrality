//Calculate Betweenness Centrality
#include <igraph.h>
#include <math.h>
#include <sys/time.h>

/**
    Get current time in milisecond.

    @return the current time in milisecond
*/
double get_current_time_mlsec() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 1e-6;
}

/**
    Print igraph_vector_t
*/
void print_vector(igraph_vector_t *v) {
    long int i;
    for (i=0; i<igraph_vector_size(v); i++) {
        fprintf(stdout, " %lf\n", VECTOR(*v)[i]);
    }
    //fprintf(f, "\n");
}


//Main function
int main() {
    igraph_t g;
    igraph_vector_t bet;
    FILE *f;

    f=fopen("../data/test2.txt", "r");
    igraph_read_graph_edgelist(&g, f, 0, 1);

    igraph_vector_init(&bet, 0);

    // igraph_betweenness(/*graph=*/ &g, /*res=*/ &bet, /*vids=*/ igraph_vss_all(),
    //          /*directed=*/1, /*weights=*/ NULL, /*nobigint=*/ 0);

    igraph_articulation_points(&g, &bet);
    igraph_vector_sort(&bet);
    print_vector(&bet);

    fclose(f);
    //graph_vector_destroy(&bet);
    igraph_destroy(&g);
    return 0;
}
