/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file libcrbbasic/mprintf.c
 *
 *  @brief Printf wrapper to allow us to control output
 *
 *  Module:  mprintf
 *
 *  Library: crbbasic
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-01-03
 *
 *
 *  Revision History:
 *         - 2008Jan03 bienert: created
 *
 */

#include <config.h>
#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
#include <libcrbfallback/snprintf.h>
#include "inc_strg.h"
#include "mprintf_dummy.h"


/** @brief Write to given output stream using a @c va_list.
 *
 * Writes a format string to a given output stream. The format string is filled
 * from the @c va_list parameter. This function does not call the va_end macro.
 * Therefore it should be called manually after each call.
 * Returns the number of characters printed (without the trailing '\0').
 * @param[in]  stream Output stream to write to.
 * @param[in]  format Format string.
 * @param[in]  ap @c va_list.
*/
int
#ifdef __GNUC__
__attribute__ ((format (printf, 2, 0)))
#endif /* __GNUC__ */
mvfprintf (FILE *stream, const char *format, va_list ap)
{
   int retval;
   
   retval = vfprintf (stream, format, ap);

   return retval;
}


/** @brief Write to given output stream.
 *
 * Writes a format string to a given output stream.
 * Returns the number of characters printed (without the trailing '\0').
 * @param[in]  stream Output stream to write to.
 * @param[in]  format Format string.
 * @param[in]  ... Argument list for the format string.
*/
int
#ifdef __GNUC__
__attribute__ ((format (printf, 2, 3)))
#endif /* __GNUC__ */
mfprintf (FILE *stream, const char *format, ...)
{
   int retval;
   va_list ap;

   va_start (ap, format);
   retval = mvfprintf (stream, format, ap);
   va_end (ap);

   return retval;
}


/** @brief Write to stdout using a @c va_list.
 *
 * Writes a format string to stdout. The format string is filled
 * from the @c va_list parameter. This function does not call the va_end macro.
 * Therefore it should be called manually after each call.
 * Returns the number of characters printed (without the trailing '\0').
 * @param[in]  format Format string.
 * @param[in]  ap @c va_list.
*/
int
#ifdef __GNUC__
__attribute__ ((format (printf, 1, 0)))
#endif /* __GNUC__ */
mvprintf (const char *format, va_list ap)
{
   int retval;

   retval = mvfprintf (stdout, format, ap);

   return retval;
}


/** @brief Write to stdout.
 *
 * Writes a format string to stdout.
 * Returns the number of characters printed (without the trailing '\0').
 * @param[in]  format Format string.
 * @param[in]  ... Argument list for the format string.
*/
int
#ifdef __GNUC__
__attribute__ ((format (printf, 1, 2)))
#endif /* __GNUC__ */
mprintf (const char *format, ...)
{
   int retval;
   va_list ap;

   va_start (ap, format);
   retval = mvprintf (format, ap);
   va_end (ap);

   return retval;
}

/** @brief Write to a string using a @c va_list.
 *
 * Writes a format string to a string. The format string is filled
 * from the @c va_list parameter. This function does not call the va_end macro.
 * Therefore it should be called manually after each call.
 * Returns the number of characters printed (without the trailing '\0').
 * @param[out] str Output string.
 * @param[in]  format Format string.
 * @param[in]  ap @c va_list.
*/
int
#ifdef __GNUC__
__attribute__ ((format (printf, 2, 0)))
#endif /* __GNUC__ */
mvsprintf (char *str, const char *format, va_list ap)
{
   int retval;

   retval = vsprintf (str, format, ap);

   return retval;
}

/** @brief Write to a string.
 *
 * Writes a format string to a string.
 * Returns the number of characters printed (without the trailing '\0').
 * @param[out] str Output string.
 * @param[in]  format Format string.
 * @param[in]  ... Argument list for the format string.
*/
int
#ifdef __GNUC__
__attribute__ ((format (printf, 2, 3)))
#endif /* __GNUC__ */
msprintf (char *str, const char *format, ...)
{
   int retval;
   va_list ap;
   
   va_start (ap, format);
   retval = mvsprintf (str, format, ap);
   va_end (ap);

   return retval;  
}


/** @brief Write to a string with length limitation using a @c va_list.
 *
 * Writes a format string to a string. The format string is filled
 * from the @c va_list parameter. This function does not call the va_end macro.
 * Therefore it should be called manually after each call.\n
 * @c mvsnprintf() does not write more than @c size bytes to @c str (including
 * the trailing '\0'). If the output was truncated (expanded format string
 * larger than @c size bytes) then the return value is the number of characters
 * (not including the trailing '\0') of a fully expanded format string.
 * Therefore, a return value of @c size or larger signals that the output was
 * truncated.\n
 * Returns the number of characters in the completely expanded format string.
 * @param[out] str Output string.
 * @param[in]  size Number of bytes to write.
 * @param[in]  format Format string.
 * @param[in]  ap @c va_list.
*/
int
#ifdef __GNUC__
__attribute__ ((format (printf, 3, 0)))
#endif /* __GNUC__ */
mvsnprintf (char *str, const size_t size, const char *format, va_list ap)
{
   int retval;

   retval = vsnprintf (str, size, format, ap);

   return retval;
}


/** @brief Write to a string with length limitation .
 *
 * Writes a format string to a string. @c msnprintf() does not write more than
 * @c size bytes to @c str (including the trailing '\0'). If the output was
 * truncated (expanded format string larger than @c size bytes) then the return
 * value is the number of characters (not including the trailing '\0') of a
 * fully expanded format string. Therefore, a return value of @c size or larger
 * signals that the output was truncated.\n
 * Returns the number of characters in the completely expanded format string.
 * @param[out] str Output string.
 * @param[in]  size Number of bytes to write.
 * @param[in]  format Format string.
 * @param[in]  ... Argument list for the format string.
*/
int
#ifdef __GNUC__
__attribute__ ((format (printf, 3, 4)))
#endif /* __GNUC__ */
msnprintf (char *str, const size_t size, const char *format, ...)
{
   int retval;
   va_list ap;
   
   va_start (ap, format);
   retval = mvsnprintf (str, size, format, ap);
   va_end (ap);   

   return retval;
}
