/*
 * Copyright (C) 2008 Stefan Bienert
 *
 * This file is part of CoRB.
 *
 * CoRB is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CoRB is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CoRB.  If not, see <http://www.gnu.org/licenses/>.
 */


/*
 ****   Documentation header   ***
 *
 *  @file libcrbbasic/errormsg.c
 *
 *  @brief Error messaging and related functions.
 *
 *  Module: errormsg
 *
 *  Library: crbbasic
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-01-08
 *
 *
 *  Revision History:
 *         - 2008Jan08 bienert: created
 *
 */

#include <config.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>
#include <errno.h>
#include <assert.h>
#include "inc_bool.h"
#include "inc_pthr.h"
#include "inc_strg.h"
#include "mprintf.h"
#include "memmgr.h"
#include "errormsg.h"

static int
throw_error_msg_internal (const char*, int, const char*, va_list)
#ifdef __GNUC__
__attribute__ ((format (printf, 3, 0)))
#endif /* __GNUC__ */
;

static int
throw_warn_msg_internal (const char*, int, const char*, va_list)
#ifdef __GNUC__
__attribute__ ((format (printf, 3, 0)))
#endif /* __GNUC__ */
;

/** @brief Takes all components for error messaging
 *
 */
typedef struct
{
      /**< Error messenger */
      int 
#if defined(__GNUC__) && (__GNUC__ > 2)
      __attribute__ ((format (printf, 3, 0)))
#endif /* __GNUC__ */
      (*error_msgr)(const char*, int, const char*, va_list);
      /**< Warning messenger */
      int
#if defined(__GNUC__) && (__GNUC__ > 2)
      __attribute__ ((format (printf, 3, 0)))
#endif /* __GNUC__ */
      (*warn_msgr)(const char*, int, const char*, va_list);
      /**< Program name used for the messages */
      char* prog_name;
#ifdef HAVE_PTHREAD
      pthread_mutex_t err_mutex;  /**< thread safety */
      pthread_mutex_t pn_mutex;  /**< thread safety */      
#endif
} ErrorInfo;


/* All information needed by the messaging functions */
static ErrorInfo __error_info = {throw_error_msg_internal, /* error_msgr */
                                 throw_warn_msg_internal, /* warn_msgr */
                                 NULL,  /* prog_name */
#ifdef HAVE_PTHREAD
                                 PTHREAD_MUTEX_INITIALIZER,
                                 PTHREAD_MUTEX_INITIALIZER
#endif
};


/** @brief Store the name of the program.
 *
 * Stores a string used by the messaging functions of the errormsg module as
 * program name. @c prog_name is stored as copy. @c prog_name may not be
 * provided as @c NULL. @c NULL is the default value of the program name, if
 * you want to reset it, use @c free_progname().\n
 * Returns 0 on success. If you are trying to set the program name, while it is
 * already set (non @c NULL), @c ERR_PNAME_SET is returned. If @c prog_name is
 * NULL or has the null byte ('\0') in the first position, @c ERR_PNAME_EMPTY
 * is returned.
 * @param[in] prog_name Input string.
 */
int
set_progname (const char* prog_name)
{
   size_t len;

   if ((prog_name == NULL)||(prog_name[0] == '\0'))
   {
      THROW_ERROR_MSG ("Attempt to set empty string as program name.");
      return ERR_NAME_EMPTY;
   }

#ifdef HAVE_PTHREAD
   pthread_mutex_lock (&__error_info.pn_mutex);
#endif

   if (__error_info.prog_name != NULL)
   {
      THROW_ERROR_MSG ("Attempt to overwrite program name.");
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock (&__error_info.pn_mutex);
#endif
      return ERR_PNAME_SET;
   }

   len = strlen (prog_name) + 1; /* + 1 for trailing '\0' */
   __error_info.prog_name = XMALLOC (len * sizeof (*__error_info.prog_name));
   strncpy (__error_info.prog_name, prog_name, len);

#ifdef HAVE_PTHREAD
   pthread_mutex_unlock (&__error_info.pn_mutex);
#endif
   return 0;
}


/** @brief Get the program name.
 *
 * Returns a pointer to the program name (NULL indicates an empty program name).
 */
char*
get_progname (void)
{
   return __error_info.prog_name;
}


/** @brief Add name to the program name.
 *
 * Extend the current program name by @c string. @c string may not be @c NULL.\n
 * Returns 0 on success, @c ERR_NAME_EMPTY if @c string is NULL or starts with
 * the null byte ('\0'). @c ERR_PNAME_EMPTY is the error code, if the program
 * name was not set before. If the length of the program name and the string to
 * be added are to long, @c ERR_PS_TO_LONG is returned.
 * @param[in] string Input string.
 */
int
add_2_progname (const char* string)
{
   size_t string_len;
   size_t pname_len;

   if((string == NULL)||(string[0] == '\0'))
   {
      THROW_ERROR_MSG ("Attempt to add empty string to program name.");
      return ERR_NAME_EMPTY;
   }

#ifdef HAVE_PTHREAD
   pthread_mutex_lock (&__error_info.pn_mutex);
#endif

   if (__error_info.prog_name == NULL)
   {
      THROW_ERROR_MSG ("Attempt to add \"%s\" to empty program name.", string);
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock (&__error_info.pn_mutex);
#endif
      return ERR_PNAME_EMPTY;
   }

   /* realloc */
   string_len = strlen (string);
   pname_len  = strlen (__error_info.prog_name);
   /* check if __error_info.program_name + string is to large */
   if (  (SIZE_MAX - (string_len * sizeof (*string))) 
       < (pname_len * sizeof (*__error_info.prog_name)))
   {
      THROW_ERROR_MSG ("Program name \"%s\" and string \"%s\" to be added are "
                       "to long to store.", get_progname(), string);

#ifdef HAVE_PTHREAD
      pthread_mutex_unlock (&__error_info.pn_mutex);
#endif
      return ERR_PS_TO_LONG;
   }
   __error_info.prog_name = XREALLOC (__error_info.prog_name,
                                     (pname_len+string_len+1)
                                      * sizeof (*__error_info.prog_name));

   /* concatenate */
   strncat (__error_info.prog_name, string, string_len);
   
#ifdef HAVE_PTHREAD
   pthread_mutex_unlock (&__error_info.pn_mutex);
#endif
   return 0;
}


/** @brief Adds a tool name to the program name.
 *
 * Extend the current program name by @c tool_name delimited by a dash
 * (' ')).\n
 * Returns the same values as @c add_2_progname().
 * @param[in] tool_name Input string.
 */
int
add_name_2_progname (const char* tool_name)
{
   int retval;
   char* tmp;

   if((tool_name == NULL)||(tool_name[0] == '\0'))
   {
      THROW_ERROR_MSG ("Attempt to add empty string to program name.");
      return ERR_NAME_EMPTY;
   }

   tmp = XMALLOC (sizeof (*tmp) * (strlen (tool_name) + 2));
   tmp[0] = ' ';
   tmp[1] = '\0';
   strcat (tmp, tool_name);

   retval = add_2_progname (tmp);

   XFREE (tmp);

   return retval;
}


/** @brief Get the length of the program name.
 *
 * Returns the length of the program name (without the trailing '\0'), 0 if
 * @c NULL.
 */
size_t
get_progname_len (void)
{
   size_t progname_len = 0;
#ifdef HAVE_PTHREAD
   pthread_mutex_lock (&__error_info.pn_mutex);
#endif

   if (__error_info.prog_name != NULL)
   {
      progname_len = strlen (__error_info.prog_name);
   }

#ifdef HAVE_PTHREAD
   pthread_mutex_unlock (&__error_info.pn_mutex);
#endif
   return progname_len;
}


/** @brief Frees the memory of program name.
 *
 * Frees the memory allocated for the program name and sets it to @c NULL.
 */
void
free_progname (void)
{
#ifdef HAVE_PTHREAD
   pthread_mutex_lock (&__error_info.pn_mutex);
   pthread_mutex_lock (&__error_info.err_mutex);
#endif

   XFREE (__error_info.prog_name);
   __error_info.prog_name = NULL;

#ifdef HAVE_PTHREAD
   pthread_mutex_unlock (&__error_info.pn_mutex);
   pthread_mutex_unlock (&__error_info.err_mutex);
   pthread_mutex_destroy(&__error_info.pn_mutex);
   pthread_mutex_destroy(&__error_info.err_mutex);
#endif
}


/** @brief Print a message primer
 *
 * Prints program name and message type to stderr (if not @c NULL).\n
 * Returns the number of characters written or a negative value.
 *
 * @param[in] type string.
 */
# define CHK_RETVAL                           \
      if ((retval < 0)||(retval_ok < 0))\
      {\
         retval_ok = -1;\
      }\
      else\
      {\
         retval_ok += retval;\
      }

static int
print_fmt_msg_primer (const char* type)
{
   char* prog_name;
   int retval;
   int retval_ok;

   retval_ok = 0;

#ifdef HAVE_PTHREAD
   if (pthread_mutex_trylock (&__error_info.pn_mutex) == 0)
   {
#endif

      prog_name = get_progname();
      if (prog_name != NULL)
      {
         retval = mfprintf (stderr, "%s", prog_name);
         CHK_RETVAL;
         
         if (type != NULL)
         {
            retval = mfprintf (stderr,":");
            CHK_RETVAL;        
         }
      }
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock (&__error_info.pn_mutex);
  }
#endif

   if (type != NULL)
   {
      retval = mfprintf (stderr,"%s", type);
      CHK_RETVAL;
   }

   return retval_ok;
}


/** @brief Print a message to stderr.
 *
 * This function writes a message of a certain format to @c stderr. If
 * @c prog_name is not @c NULL, the output starts with "@c prog_name:". If
 * @c type is not @c NULL, the message proceeds with "@c type:". If @c file is
 * given, @c type is followed by ":@c file:@c line:". The @c format string is
 * evaluated to the message and put in the next position of the output. If the
 * message ends with a colon ":" or @c format is @c NULL, the output string is
 * extended by "@c strerror(errno)". Writing a message without @c format and
 * @c strerror(errno) can be realised by using ("%c", '\0') as format string
 * and parameter list. The string written is followed by a single newline.\n
 * Returns the number of characters printed (without the trailing '\0') or a
 * negative value on failure.
 * @param[in]  type @c type string for the message.
 * @param[in]  file Calling file name.
 * @param[in]  line Calling line.
 * @param[in]  format Format string.
 * @param[in]  ap Argument list for the format string.
 */
static int
#ifdef __GNUC__
__attribute__ ((format (printf, 4, 0)))
#endif /* __GNUC__ */
print_fmt_msg (const char* type,
               const char* file, const int line,
               const char* format,
               va_list ap)
{
   int fflushret;
   int retval = 0;
   int retval_ok = 0;

   /* clean stdout */
   fflushret = fflush (stdout);

   if (fflushret == EOF)
   {
      retval = print_fmt_msg_primer(type);
      CHK_RETVAL;

      retval = mfprintf (stderr, "Problem at writing message: %s\n",
                          strerror (errno));
      CHK_RETVAL;

      retval = mfprintf (stderr, "This is the original message:\n");
      CHK_RETVAL;
   }

   /* print prog_name:type */
   retval = print_fmt_msg_primer(type);   
   CHK_RETVAL;

   /* print file:line */
   if (file != NULL)
   {
      if (retval > 0)
      {
         retval = mfprintf (stderr,":");
         CHK_RETVAL; 
      }

      retval = mfprintf (stderr, "%s:%d", file, line);
      CHK_RETVAL;
   }

   /* print message */
   if ((format != NULL)&&(format[0] != '\0'))
   {
      if (retval > 0)
      {
         retval = mfprintf (stderr,":");
         CHK_RETVAL; 
      }
      retval = mvfprintf (stderr, format, ap);
      CHK_RETVAL;
   }

   /* print message from glibc, if desired */
   if ((format == NULL)
       ||((format[0] != '\0')&&(format[strlen(format) - 1] == ':')))
   {
      retval = mfprintf (stderr, "%s", strerror(errno));
      CHK_RETVAL;
   }
   retval = mfprintf (stderr, "\n");
   CHK_RETVAL;

   return retval_ok;
}


/** @brief Print an error message to @c stderr.
 *
 * Writes a message with the word "ERROR" in it to @c stderr. If the program
 * name is set, it is used as a primer to the message. If @c file is
 * provided, @c file and @line are added to the message. If @c format ends on a
 * colon or is @c NULL, @c strerror(errno) is added to the end of the message.
 * No message at all is produced, if ("%c", '\0') is used as format string and
 * parameter list.\n
 * Returns the return value of @c print_fmt_msg().\n
 * If used via the macro @c THROW_ERROR_MSG, @c file and @c line are set
 * automatically when compiled without @c NDEBUG.
 * @param[in] file Calling file. Usually filled by the __FILE__ macro.
 * @param[in] line Called from line. Use __LINE__ here.
 * @param[in] format Format string.
 * @param[in] ap   Parameter list.
 */
static int
throw_error_msg_internal (const char* file, const int line,
                          const char* format,
                          va_list ap)
{
   return print_fmt_msg ("ERROR", file, line, format, ap);
}


/** @brief Get the current error message function
 *
 * This function should return a pointer to the message function. For the case
 * that the functionw as set to @c NULL, @c NULL is returned WITHOUT an error
 * message because of lack of an error message function.\n
 * Returns a pointer to the message function.
 */
int
(*get_error_msg_func(void))(const char*, int, const char*, va_list)
{
   int ((*msg_func)(const char*, int, const char*, va_list));

#ifdef HAVE_PTHREAD
   pthread_mutex_lock (&__error_info.err_mutex);
#endif

   assert (__error_info.error_msgr != NULL);

   msg_func = __error_info.error_msgr;

#ifdef HAVE_PTHREAD
   pthread_mutex_unlock (&__error_info.err_mutex);
#endif
   return msg_func;
}


/** @brief Set error message function
 *
 * The function to be set as error message function has to be of the following
 * prototype:\n
 * @c int (*error_msgr)(const char*, int, const char*, va_list)
 * @param[in] error_msgr Error message function to use.
 */
void
set_error_msg_func (int (*error_msgr)(const char*, int, const char*, va_list))
{
#ifdef HAVE_PTHREAD
   pthread_mutex_lock (&__error_info.err_mutex);
#endif
   __error_info.error_msgr = error_msgr;
#ifdef HAVE_PTHREAD
   pthread_mutex_unlock (&__error_info.err_mutex);
#endif
}


/** @brief Call the error message function
 *
 * This function presents one way to call the error message function. It might
 * take a file name and line number (not displayed if file is @c NULL). Your
 * own message comes with the format string and its parameters. If the format
 * string ends on a colon, actually on ":\0", or is @c NULL, then the output of
 * @c strerror(errno) is printed to @c stderr. No message string is produced,
 * if ("%c", '\0') is used as format string and parameter list. If a program
 * name is set (@c set_progname()), it will be printed at the very beginning of
 * the message.\n
 * Returns a negative value on failure, the number of written characters,
 * else.\n
 * An other way to access the error message function is via the macro
 * @c THROW_ERROR_MSG. This handles the @c file and @c line parameter for you:
 * If compiled without NDEBUG, the calling file and line will be included in
 * the message.
 * @param[in] file Calling file. Usually filled by the __FILE__ macro.
 * @param[in] line Called from line. Use __LINE__ here.
 * @param[in] format Format string.
 * @param[in] ...   Parameter list for the format string.
 */
int
call_error_msgr (const char* file, int line, const char* format, ...)
{
   int retval;
   va_list ap;

   va_start (ap, format);
   retval = get_error_msg_func() (file, line, format, ap);
   va_end(ap);

   return retval;
}

/* warning messages */
/** @brief Print an warning message to @c stderr.
 *
 * Writes a message with the word "Warning" in it to @c stderr. If the program
 * name is set, it is used as a primer to the message. If @c file is
 * provided, @c file and @line are added to the message. If @c format ends on a
 * colon or is @c NULL, @c strerror(errno) is added to the end of the message.
 * No message at all is produced, if ("%c", '\0') is used as format string and
 * parameter list.\n
 * Returns the return value of @c print_fmt_msg().\n
 * If used via the macro @c THROW_ERROR_MSG, @c file and @c line are set
 * automatically when compiled without @c NDEBUG.
 * @param[in] file Calling file. Usually filled by the __FILE__ macro.
 * @param[in] line Called from line. Use __LINE__ here.
 * @param[in] format Format string.
 * @param[in] ap   Parameter list.
 */
static int
throw_warn_msg_internal (const char* file, const int line,
                         const char* format,
                         va_list ap)
{
   return print_fmt_msg ("WARNING", file, line, format, ap);
}


/** @brief Get the current warning message function
 *
 * This function should return a pointer to the message function. For the case
 * that the functionw as set to @c NULL, @c NULL is returned WITHOUT an error
 * message.\n
 * Returns a pointer to the message function.
 */
int
(*get_warn_msg_func(void))(const char*, int, const char*, va_list)
{
   int (*msg_func)(const char*, int, const char*, va_list);

#ifdef HAVE_PTHREAD
   pthread_mutex_lock (&__error_info.err_mutex);
#endif

   assert (__error_info.warn_msgr != NULL);

   msg_func = __error_info.warn_msgr;

#ifdef HAVE_PTHREAD
   pthread_mutex_unlock (&__error_info.err_mutex);
#endif

   return msg_func;
}


/** @brief Set warning message function
 *
 * The function to be set as warning message function has to be of the following
 * prototype:\n
 * @c int (*warn_msgr)(const char*, int, const char*, va_list)
 * @param[in] warn_msgr Error message function to use.
 */
void
set_warn_msg_func (int (*warn_msgr)(const char*, int, const char*, va_list))
{
#ifdef HAVE_PTHREAD
   pthread_mutex_lock (&__error_info.err_mutex);
#endif
   __error_info.warn_msgr = warn_msgr;
#ifdef HAVE_PTHREAD
   pthread_mutex_unlock (&__error_info.err_mutex);
#endif
}


/** @brief Call the warning message function
 *
 * This function presents one way to call the warning message function. It might
 * take a file name and line number (not displayed if file is @c NULL). Your
 * own message comes with the format string and its parameters. If the format
 * string ends on a colon, actually on ":\0", or is @c NULL, then the output of
 * @c strerror(errno) is printed to @c stderr. No message string is produced,
 * if ("%c", '\0') is used as format string and parameter list. If a program
 * name is set (@c set_progname()), it will be printed at the very beginning of
 * the message.\n
 * Returns a negative value on failure, the number of written characters,
 * else.\n
 * An other way to access the warning message function is via the macro
 * @c THROW_WARN_MSG. This handles the @c file and @c line parameter for you:
 * If compiled without NDEBUG, the calling file and line will be included in
 * the message.
 * @param[in] file Calling file. Usually filled by the __FILE__ macro.
 * @param[in] line Called from line. Use __LINE__ here.
 * @param[in] format Format string.
 * @param[in] ...   Parameter list for the format string.
 */
int
call_warn_msgr (const char* file, int line, const char* format, ...)
{
   int retval;
   va_list ap;

   va_start (ap, format);
   retval = get_warn_msg_func() (file, line, format, ap);
   va_end(ap);

   return retval;
}
