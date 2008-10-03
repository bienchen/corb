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


#undef realloc

#include <config.h>
#include <sys/types.h>

void *realloc ();

void *
rpl_realloc (void* ptr, size_t n)
{
   if (n == 0)
   {
      n = 1;
      
      free (ptr);
      ptr = NULL;
   }

   return (ptr == NULL ? malloc (n) : realloc (ptr, n));
}
