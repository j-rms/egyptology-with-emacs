(defun insert-egyptian-hieroglyph ()
  "Use Helm to select and insert a Unicode Egyptian Hieroglyph at point, showing Gardiner sign number and description."
  (interactive)
  (require 'helm)
  (let* ((hieroglyphs
          (cl-loop for code from #x13000 to #x1342E
                   for name = (get-char-code-property code 'name)
                   when name
                   collect (cons (format "%c %s" code name code) code)))
         (selected (helm :sources
                         `((name . "Egyptian Hieroglyphs")
                           (candidates . ,hieroglyphs)
                           (action . (lambda (code) code)))))
         (char (when selected (char-to-string selected))))
    (when char
      (insert char))))

(global-set-key (kbd "C-c g") 'insert-egyptian-hieroglyph)
