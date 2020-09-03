## Anpassung der subop_2 Routinen

* offene Frage
  * Wirksamkeit der SIMD Direktive bei variabler Schleifengröße
  * SIMD Direktive ggf. von `i`-  auf `j`-Schleife umsetzen

### Ausgangsdateien

      ../../tpo_aaa/subop_2__llll.f
      ../../tpo_aaa/subop_2__lu2ll_r.f
      ../../tpo_aaa/subop_2__lu4ll_r.f
      ../../tpo_aaa/subop_2__lu8ll_r.f
      ../../tpo_aaa/subop_2__lu4lb4_r.f 

### Aktionen

__noch anpassen__

* Header anpassen

* Intent von `v` zu `inout` ändern

* Dimension von `u` und `v` setzen

    - erste auf `na`
    - dritte auf  `nc`
    
* Deklaration von Input-Argumenten einfügen

      integer,   intent(in)    :: nb                   !< 2nd dimension of u,v
      integer,   intent(in)    :: nc                   !< 3rd dimension of u,v
      real(RNP), intent(in)    :: alpha                !< factor α
      real(RNP), intent(in)    :: beta                 !< factor β

* Ersetzungen

      __NA1__  →  _NB2_
      __NA2__  →  _NB1_
      At       →  B
      z1       →  u
      z2       →  v
      subroutine SubOp_2(At, z1, z2)
        → subroutine PROC(IxBtxI__,_NB1_,_NB2_)(na, nc, B, alpha, beta, u, v)
      end subroutine SubOp_2
        →  end subroutine PROC(IxBtxI__,_NB1_,_NB2_)
  etc.

* Zuweisung(en) anpassen
  
      z2(i,j,k) = tmp  →  v(i,j,k) = alpha * tmp + beta * v(i,j,k)
  
  o.Ä.

### Testprogram

0. Zu testende Routine durch Anpassen der `#include` Direktive in `vec_test.F` einfügen

1. Preprozessor-Ausgabe prüfen

        ifort -free -E -P vec_test.fpp > vec_test.f

2. Korrektheit prüfen 

        ifort -free -O2 -check all -o vec_test vec_test.fpp
        ./vec_test

3. Vektorisierung anschauen

        ifort -free -O2 -o vec_test -qopt-report=2 -qopt-report-phase=vec vec_test.fpp

   * Optimierungsreport `vec_test.optrpt` im Editor öffnen
   * Mit `vec_test.fpp` sowie `vec_test.f` vergleichen
   * Hinweisen nachgehen

4. Leistungstest

        ifort -free -O2 -o vec_test vec_test.fpp
        ./vec_test

5. Änderungen vornehmen und mit 1.  fortsetzen



