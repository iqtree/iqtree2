
/**
 interface to call booster for transfer bootstrap expectation (TBE)
 @param input_tree reference tree file
 @param boot_trees bootstrap trees file
 @param out_tree output tree
 @param out_raw_tree output raw tree
 @param stat_out statistic output file
 @param num_threads number of threads
 @param quiet 1 to stay quiet, 0 otherwise
 */
int main_booster (const char* input_tree, const char *boot_trees,
                  const char* out_tree, const char* out_raw_tree, const char* stat_out,
                  int quiet);
