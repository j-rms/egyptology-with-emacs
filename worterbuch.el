;; Functions for smoothly opening Worterbuch and CT volumes

(setf bks/worterbuch-vols '((bd1 . "/home/joel/projects/library/Worterbuch/bd1-trimmed.pdf")
                            (bd2 . "/home/joel/projects/library/Worterbuch/Bd_II.pdf")
                            (bd3 . "/home/joel/projects/library/Worterbuch/Bd_III.pdf")
                            (bd4 . "/home/joel/projects/library/Worterbuch/Bd_IV.pdf")
                            (bd5 . "/home/joel/projects/library/Worterbuch/Bd_V.pdf")))

(setf bks/worterbuch-vols-page-offsets '(()))

(defun bds/get-vol-path (liblist volnum)
  (cdr (nth (- volnum 1) liblist )))

(defun Wb (volnum)
  (interactive "nWb volnum: ")
  (let ((large-file-warning-threshold nil))
    (find-file (bds/get-vol-path bks/worterbuch-vols volnum))))

(setf bks/CT-vols '((ct1 . "/home/joel/projects/library/coffin texts/oip34 vol1.pdf")
                    (ct2 . "/home/joel/projects/library/coffin texts/oip49 vol2.pdf")
                    (ct3 . "/home/joel/projects/library/coffin texts/oip64 vol3.pdf")
                    (ct4 . "/home/joel/projects/library/coffin texts/oip67 vol4.pdf")
                    (ct5 . "/home/joel/projects/library/coffin texts/oip73 vol5.pdf")
                    (ct6 . "/home/joel/projects/library/coffin texts/oip81 vol6.pdf")
                    (ct7 . "/home/joel/projects/library/coffin texts/oip87 vol7.pdf")))

(defun CT (volnum)
  (interactive "nCT volnum: ")
  (let ((large-file-warning-threshold nil))
    (find-file (bds/get-vol-path bks/CT-vols volnum))))


(setf bks/CT-faulkner-vols '((ctf1 . "/home/joel/projects/library/faulkner2004a.pdf")
                             (ctf2 . "/home/joel/projects/library/faulkner2004b.pdf")
                             (ctf3 . "/home/joel/projects/library/faulkner2004c.pdf")))

(defun CTf (volnum)
  (interactive "nCT Faulkner translation volnum: ")
  (let ((large-file-warning-threshold nil))
    (find-file (bds/get-vol-path bks/CT-faulkner-vols volnum))))
