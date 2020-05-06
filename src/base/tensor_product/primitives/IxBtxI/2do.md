## Anpassung der subop_2 Routinen

* diese Datei ab Aktionen überarbeiten
* offene Frage
  * Wirksamkeit der SIMD Direktive bei variabler Schleifengröße
  * SIMD Direktive ggf. von `i`-  auf `j`-Schleife umsetzen

### Ausgangsdateien

      ../tpo_aaa/subop_2__llll.f
      ../tpo_aaa/subop_2__lu2ll_r.f
      ../tpo_aaa/subop_2__lu4ll_r.f
      ../tpo_aaa/subop_2__lu8ll_r.f
      ../tpo_aaa/subop_2__lu4lb4_r.f 

### Aktionen

__noch anpassen__

* Header anpassen

* Dimension von `u` und `v` setzen

    - erste auf `na`
    - dritte auf  `nc`
    
* Deklaration von Input-Argumenten einfügen

      integer,   intent(in)    :: nb                   !< 2nd dimension of u,v
      integer,   intent(in)    :: nc                   !< 3rd dimension of u,v
      real(RNP), intent(in)    :: alpha                !< factor α
      real(RNP), intent(in)    :: beta                 !< factor β

* Ersetzungen

      __NA1__  →  _NA2_
      __NA2__  →  _NA1_
      At       →  A
      z1       →  v
      subroutine SubOp_1(At, u, z1)
        → subroutine PROC(IxIxAt__,_NA1_,_NA2_)(A, nb, nc, alpha, beta, u, v)
      end subroutine SubOp_1
        →  end subroutine PROC(IxIxAt__,_NA1_,_NA2_)

* Zuweisung(en) anpassen
  
      z1(i,j,k) = tmp  →  v(i,j,k) = alpha * tmp + beta * v(i,j,k)
  

o.ä.
