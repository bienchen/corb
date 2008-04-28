;; emacs font-lock enhancements for libgtcore

;; Array, Str, Env, GenFile
(font-lock-add-keywords
 'c-mode
 '(
   ("\\<\\(Array\\)\\>" 1 font-lock-type-face)
   ("\\<\\(Str\\)\\>" 1 font-lock-type-face)
   ("\\<\\(Env\\)\\>" 1 font-lock-type-face)
   ("\\<\\(GenFile\\)\\>" 1 font-lock-type-face)
   ))
