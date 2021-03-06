/** @file brot_cmdline.h
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.22.1
 *  http://www.gnu.org/software/gengetopt.
 *  DO NOT modify this file, since it can be overwritten
 *  @author GNU Gengetopt by Lorenzo Bettini */

#ifndef BROT_CMDLINE_H
#define BROT_CMDLINE_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef BROT_CMDLINE_PARSER_PACKAGE
/** @brief the program name */
#define BROT_CMDLINE_PARSER_PACKAGE "brot"
#endif

#ifndef BROT_CMDLINE_PARSER_VERSION
/** @brief the program version */
#define BROT_CMDLINE_PARSER_VERSION "1.0"
#endif

enum enum_scoring { scoring_arg_NN = 0 , scoring_arg_nussinov, scoring_arg_simpleNN };

/** @brief Where the command line options are stored */
struct brot_args_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *detailed_help_help; /**< @brief Print help, including all details and hidden options, and exit help description.  */
  const char *full_help_help; /**< @brief Print help, including hidden options, and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  const char *verbose_help; /**< @brief Enable verbose mode help description.  */
  const char *file_help; /**< @brief Read structure from file help description.  */
  enum enum_scoring scoring_arg;	/**< @brief Choose a different scoring scheme (default='NN').  */
  char * scoring_orig;	/**< @brief Choose a different scoring scheme original value given at command line.  */
  const char *scoring_help; /**< @brief Choose a different scoring scheme help description.  */
  char ** fixed_nuc_arg;	/**< @brief Preset a nucleotide in a position.  */
  char ** fixed_nuc_orig;	/**< @brief Preset a nucleotide in a position original value given at command line.  */
  unsigned int fixed_nuc_min; /**< @brief Preset a nucleotide in a position's minimum occurreces */
  unsigned int fixed_nuc_max; /**< @brief Preset a nucleotide in a position's maximum occurreces */
  const char *fixed_nuc_help; /**< @brief Preset a nucleotide in a position help description.  */
  long steps_arg;	/**< @brief Number of iterations (default='1000').  */
  char * steps_orig;	/**< @brief Number of iterations original value given at command line.  */
  const char *steps_help; /**< @brief Number of iterations help description.  */
  float temp_arg;	/**< @brief Initial temperature (default='2').  */
  char * temp_orig;	/**< @brief Initial temperature original value given at command line.  */
  const char *temp_help; /**< @brief Initial temperature help description.  */
  long seed_arg;	/**< @brief Random seed.  */
  char * seed_orig;	/**< @brief Random seed original value given at command line.  */
  const char *seed_help; /**< @brief Random seed help description.  */
  float negative_design_scaling_arg;	/**< @brief Scale negative design term (default='0.42').  */
  char * negative_design_scaling_orig;	/**< @brief Scale negative design term original value given at command line.  */
  const char *negative_design_scaling_help; /**< @brief Scale negative design term help description.  */
  float heterogenity_term_scaling_arg;	/**< @brief Scale heterogenity term (default='9.73').  */
  char * heterogenity_term_scaling_orig;	/**< @brief Scale heterogenity term original value given at command line.  */
  const char *heterogenity_term_scaling_help; /**< @brief Scale heterogenity term help description.  */
  char * entropy_output_arg;	/**< @brief Write down entropy and temperature changes.  */
  char * entropy_output_orig;	/**< @brief Write down entropy and temperature changes original value given at command line.  */
  const char *entropy_output_help; /**< @brief Write down entropy and temperature changes help description.  */
  char * simulation_output_arg;	/**< @brief Write down the matrix of each step.  */
  char * simulation_output_orig;	/**< @brief Write down the matrix of each step original value given at command line.  */
  const char *simulation_output_help; /**< @brief Write down the matrix of each step help description.  */
  long window_size_arg;	/**< @brief Window size for the heterogeneity term (default='1').  */
  char * window_size_orig;	/**< @brief Window size for the heterogeneity term original value given at command line.  */
  const char *window_size_help; /**< @brief Window size for the heterogeneity term help description.  */
  float sm_entropy_arg;	/**< @brief Sequence matrix entropy threshold (default='0.337').  */
  char * sm_entropy_orig;	/**< @brief Sequence matrix entropy threshold original value given at command line.  */
  const char *sm_entropy_help; /**< @brief Sequence matrix entropy threshold help description.  */
  float lambda_arg;	/**< @brief Portion of a new step to be accepted (default='0.627').  */
  char * lambda_orig;	/**< @brief Portion of a new step to be accepted original value given at command line.  */
  const char *lambda_help; /**< @brief Portion of a new step to be accepted help description.  */
  float beta_long_arg;	/**< @brief Define the long term entropy contribution (default='0.949').  */
  char * beta_long_orig;	/**< @brief Define the long term entropy contribution original value given at command line.  */
  const char *beta_long_help; /**< @brief Define the long term entropy contribution help description.  */
  float beta_short_arg;	/**< @brief Define the short term entropy contribution (default='0.5').  */
  char * beta_short_orig;	/**< @brief Define the short term entropy contribution original value given at command line.  */
  const char *beta_short_help; /**< @brief Define the short term entropy contribution help description.  */
  float speedup_threshold_arg;	/**< @brief Speedup/ slow down cooling threshold (default='0.816').  */
  char * speedup_threshold_orig;	/**< @brief Speedup/ slow down cooling threshold original value given at command line.  */
  const char *speedup_threshold_help; /**< @brief Speedup/ slow down cooling threshold help description.  */
  float min_cool_arg;	/**< @brief Minimal cooling factor (default='0.866').  */
  char * min_cool_orig;	/**< @brief Minimal cooling factor original value given at command line.  */
  const char *min_cool_help; /**< @brief Minimal cooling factor help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int detailed_help_given ;	/**< @brief Whether detailed-help was given.  */
  unsigned int full_help_given ;	/**< @brief Whether full-help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int verbose_given ;	/**< @brief Whether verbose was given.  */
  unsigned int file_given ;	/**< @brief Whether file was given.  */
  unsigned int scoring_given ;	/**< @brief Whether scoring was given.  */
  unsigned int fixed_nuc_given ;	/**< @brief Whether fixed-nuc was given.  */
  unsigned int steps_given ;	/**< @brief Whether steps was given.  */
  unsigned int temp_given ;	/**< @brief Whether temp was given.  */
  unsigned int seed_given ;	/**< @brief Whether seed was given.  */
  unsigned int negative_design_scaling_given ;	/**< @brief Whether negative-design-scaling was given.  */
  unsigned int heterogenity_term_scaling_given ;	/**< @brief Whether heterogenity-term-scaling was given.  */
  unsigned int entropy_output_given ;	/**< @brief Whether entropy-output was given.  */
  unsigned int simulation_output_given ;	/**< @brief Whether simulation-output was given.  */
  unsigned int window_size_given ;	/**< @brief Whether window-size was given.  */
  unsigned int sm_entropy_given ;	/**< @brief Whether sm-entropy was given.  */
  unsigned int lambda_given ;	/**< @brief Whether lambda was given.  */
  unsigned int beta_long_given ;	/**< @brief Whether beta-long was given.  */
  unsigned int beta_short_given ;	/**< @brief Whether beta-short was given.  */
  unsigned int speedup_threshold_given ;	/**< @brief Whether speedup-threshold was given.  */
  unsigned int min_cool_given ;	/**< @brief Whether min-cool was given.  */

  char **inputs ; /**< @brief unamed options (options without names) */
  unsigned inputs_num ; /**< @brief unamed options number */
} ;

/** @brief The additional parameters to pass to parser functions */
struct brot_cmdline_parser_params
{
  int override; /**< @brief whether to override possibly already present options (default 0) */
  int initialize; /**< @brief whether to initialize the option structure brot_args_info (default 1) */
  int check_required; /**< @brief whether to check that all required options were provided (default 1) */
  int check_ambiguity; /**< @brief whether to check for options already specified in the option structure brot_args_info (default 0) */
  int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
} ;

/** @brief the purpose string of the program */
extern const char *brot_args_info_purpose;
/** @brief the usage string of the program */
extern const char *brot_args_info_usage;
/** @brief all the lines making the help output */
extern const char *brot_args_info_help[];
/** @brief all the lines making the full help output (including hidden options) */
extern const char *brot_args_info_full_help[];
/** @brief all the lines making the detailed help output (including hidden options and details) */
extern const char *brot_args_info_detailed_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int brot_cmdline_parser (int argc, char **argv,
  struct brot_args_info *args_info);

/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use brot_cmdline_parser_ext() instead
 */
int brot_cmdline_parser2 (int argc, char **argv,
  struct brot_args_info *args_info,
  int override, int initialize, int check_required);

/**
 * The command line parser (version with additional parameters)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int brot_cmdline_parser_ext (int argc, char **argv,
  struct brot_args_info *args_info,
  struct brot_cmdline_parser_params *params);

/**
 * Save the contents of the option struct into an already open FILE stream.
 * @param outfile the stream where to dump options
 * @param args_info the option struct to dump
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int brot_cmdline_parser_dump(FILE *outfile,
  struct brot_args_info *args_info);

/**
 * Save the contents of the option struct into a (text) file.
 * This file can be read by the config file parser (if generated by gengetopt)
 * @param filename the file where to save
 * @param args_info the option struct to save
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int brot_cmdline_parser_file_save(const char *filename,
  struct brot_args_info *args_info);

/**
 * Print the help
 */
void brot_cmdline_parser_print_help(void);
/**
 * Print the full help (including hidden options)
 */
void brot_cmdline_parser_print_full_help(void);
/**
 * Print the detailed help (including hidden options and details)
 */
void brot_cmdline_parser_print_detailed_help(void);
/**
 * Print the version
 */
void brot_cmdline_parser_print_version(void);

/**
 * Initializes all the fields a brot_cmdline_parser_params structure 
 * to their default values
 * @param params the structure to initialize
 */
void brot_cmdline_parser_params_init(struct brot_cmdline_parser_params *params);

/**
 * Allocates dynamically a brot_cmdline_parser_params structure and initializes
 * all its fields to their default values
 * @return the created and initialized brot_cmdline_parser_params structure
 */
struct brot_cmdline_parser_params *brot_cmdline_parser_params_create(void);

/**
 * Initializes the passed brot_args_info structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void brot_cmdline_parser_init (struct brot_args_info *args_info);
/**
 * Deallocates the string fields of the brot_args_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void brot_cmdline_parser_free (struct brot_args_info *args_info);

/**
 * The string parser (interprets the passed string as a command line)
 * @param cmdline the command line stirng
 * @param args_info the structure where option information will be stored
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int brot_cmdline_parser_string (const char *cmdline, struct brot_args_info *args_info,
  const char *prog_name);
/**
 * The string parser (version with additional parameters - deprecated)
 * @param cmdline the command line stirng
 * @param args_info the structure where option information will be stored
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use brot_cmdline_parser_string_ext() instead
 */
int brot_cmdline_parser_string2 (const char *cmdline, struct brot_args_info *args_info,
  const char *prog_name,
  int override, int initialize, int check_required);
/**
 * The string parser (version with additional parameters)
 * @param cmdline the command line stirng
 * @param args_info the structure where option information will be stored
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int brot_cmdline_parser_string_ext (const char *cmdline, struct brot_args_info *args_info,
  const char *prog_name,
  struct brot_cmdline_parser_params *params);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int brot_cmdline_parser_required (struct brot_args_info *args_info,
  const char *prog_name);

extern const char *brot_cmdline_parser_scoring_values[];  /**< @brief Possible values for scoring. */


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* BROT_CMDLINE_H */
