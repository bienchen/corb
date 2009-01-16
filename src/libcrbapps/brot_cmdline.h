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
#define BROT_CMDLINE_PARSER_VERSION "0.1"
#endif

enum enum_scoring { scoring_arg_NN = 0 , scoring_arg_nussinov, scoring_arg_simpleNN };

/** @brief Where the command line options are stored */
struct brot_args_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *detailed_help_help; /**< @brief Print help, including all details and hidden options, and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
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
  float temp_arg;	/**< @brief Initial temperature (default='100').  */
  char * temp_orig;	/**< @brief Initial temperature original value given at command line.  */
  const char *temp_help; /**< @brief Initial temperature help description.  */
  long seed_arg;	/**< @brief Random seed (default='-791122').  */
  char * seed_orig;	/**< @brief Random seed original value given at command line.  */
  const char *seed_help; /**< @brief Random seed help description.  */
  float negative_design_scaling_arg;	/**< @brief Scale negative design term (default='1.000000').  */
  char * negative_design_scaling_orig;	/**< @brief Scale negative design term original value given at command line.  */
  const char *negative_design_scaling_help; /**< @brief Scale negative design term help description.  */
  float heterogenity_term_scaling_arg;	/**< @brief Scale heterogenity term (default='1.000000').  */
  char * heterogenity_term_scaling_orig;	/**< @brief Scale heterogenity term original value given at command line.  */
  const char *heterogenity_term_scaling_help; /**< @brief Scale heterogenity term help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int detailed_help_given ;	/**< @brief Whether detailed-help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int scoring_given ;	/**< @brief Whether scoring was given.  */
  unsigned int fixed_nuc_given ;	/**< @brief Whether fixed-nuc was given.  */
  unsigned int steps_given ;	/**< @brief Whether steps was given.  */
  unsigned int temp_given ;	/**< @brief Whether temp was given.  */
  unsigned int seed_given ;	/**< @brief Whether seed was given.  */
  unsigned int negative_design_scaling_given ;	/**< @brief Whether negative-design-scaling was given.  */
  unsigned int heterogenity_term_scaling_given ;	/**< @brief Whether heterogenity-term-scaling was given.  */

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
/** @brief all the lines making the detailed help output (including hidden options and details) */
extern const char *brot_args_info_detailed_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int brot_cmdline_parser (int argc, char** argv,
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
int brot_cmdline_parser2 (int argc, char** argv,
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
int brot_cmdline_parser_ext (int argc, char** argv,
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
