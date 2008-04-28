;; emacs font-lock enhancements for the Pdb classes

;; PdbAtom
(font-lock-add-keywords
 'c-mode
 '(
   ("\\<\\(PdbAtom\\)\\>" 1 font-lock-type-face)
   ))

;; pdb_scene
(font-lock-add-keywords
 'c-mode
 '(
   ("\\<\\(Pdb_scene\\)\\>" 1 font-lock-type-face)
   ("\\<\\(Pdb_parser_state\\)\\>" 1 font-lock-type-face)
   ))


