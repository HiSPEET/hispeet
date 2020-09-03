## Anpassung der subop_3 Routinen

### Ausgangsdateien

      ../../tpo_aaa/subop_3__lffl.f
      ../../tpo_aaa/subop_3__u2ffl_r.f
      ../../tpo_aaa/subop_2__lu4ll_r.f
      ../../tpo_aaa/subop_3__u4ffl_r.f
      ../../tpo_aaa/subop_3__u8ffl_r.f
      ../../tpo_aaa/subop_3__u4ffb4_r.f

### Aktionen

__noch anpassen__

* Header anpassen

* Erste Dimension von `u` und `v`  auf `na_nb` setzen

* Intent von `v` zu `inout` ändern

* Deklaration von Input-Argumenten einfügen

      integer,   intent(in)    :: na_nb           !< size of dimensions 1+2 of u,v
      integer,   intent(in)    :: nc              !< 3rd dimension of u,v
      real(RNP), intent(in)    :: alpha           !< factor α
      real(RNP), intent(in)    :: beta            !< factor β

* Ersetzungen

      __NA1__  →  _NC2_
      __NA2__  →  _NC1_
      At       →  C
      z1       →  u
      z2       →  v
      subroutine SubOp_2(At, z2, v)
        → subroutine PROC(CtxIxI__,_NC1_,_NC2_)(na_nb, C, alpha, beta, u, v)
      end subroutine SubOp_2
        →  end subroutine PROC(CtxIxI__,_NC1_,_NC2_)
  etc.

* Zuweisung(en) anpassen
  
      v(i,j,k) = tmp  →  v(i,j,k) = alpha * tmp + beta * v(i,j,k)
  
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



