;; emacs font-lock enhancements for libgtcore

;; Copyright (C) 2007 Stefan Bienert
;; 
;; This file is part of CoRB.
;; 
;; CoRB is free software: you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation, either version 3 of the License, or
;; (at your option) any later version.
;; 
;; CoRB is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.
;; 
;; You should have received a copy of the GNU General Public License
;; along with CoRB.  If not, see <http://www.gnu.org/licenses/>.

;; Array, Str, Env, GenFile
(font-lock-add-keywords
 'c-mode
 '(
   ("\\<\\(Array\\)\\>" 1 font-lock-type-face)
   ("\\<\\(Str\\)\\>" 1 font-lock-type-face)
   ("\\<\\(Env\\)\\>" 1 font-lock-type-face)
   ("\\<\\(GenFile\\)\\>" 1 font-lock-type-face)
   ))
