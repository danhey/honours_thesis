(define-param eps-hi 12.25)	; Waveguide permittivity
(define-param eps-lo 1)		; Surrounding permittivity
(define-param h 1.1)		; Thickness (microns)
(define-param Y 5) 			; Computational size Y

; Generate lattice
(set! geometry-lattice (make lattice (size no-size Y no-size)))
(set! default-material (make dielectric (epsilon eps-lo)))

;Implement waveguide structure
(set! geometry
 (list (make block 			; A dielectric,
 (center 0 0 0) 			; centered at the origin,
 (size infinity h infinity) ; of infinite length
 (material (make dielectric (epsilon eps-hi))))))
(set-param! resolution 50)

;k-points over which to interpolate
(define-param kmin 2)
(define-param kmax 3)
(define-param k-interp 200)

(set! k-points (interpolate k-interp
 (list (vector3 kmin 0 0) (vector3 kmax 0 0))))
 
;Number of bands to calculate
(set-param! num-bands 2)

(run-te-yeven)
(run-te-yodd)