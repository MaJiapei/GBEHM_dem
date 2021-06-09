! A fortran95 program for G95
! By WQY
module Priver
    implicit none
    Character (100)::gis_dir="./", para_dir="./"

contains
        Subroutine order(ac, bs, bn)
          implicit none
          Integer, Parameter :: maxb = 500
          Integer ac(maxb)
          Integer bs(maxb)
          Integer bn
          Integer i, j, tmp
          Do i = 1, bn - 1
            Do j = i + 1, bn
              If (ac(j)<ac(i)) Then
                tmp = bs(i)
                bs(i) = bs(j)
                bs(j) = tmp
                tmp = ac(i)
                ac(i) = ac(j)
                ac(j) = tmp
              End If
            End Do
          End Do
        End Subroutine order

        Subroutine links(dir, netbasinend, net, nsub, nbasin, basinpnext, outlet)
          Integer dir(1000, 1000), net(1000, 1000)
          Integer netbasinend(1000, 2), basinpnext(1000)
          Integer nbasin(1000)
          Integer nsub
          Integer i, j, k, m
          Integer nexti, nextj
          Integer outlet
          basinpnext = -9999
          Do k = 1, nsub
            If (k/=outlet) Then
              m = 0
              i = netbasinend(k, 1)
              j = netbasinend(k, 2)
              Do While (m==0)
                Print *, i, j
                If (dir(i,j)==64) Then
                  nexti = i - 1
                  nextj = j
                Else If (dir(i,j)==128) Then
                  nexti = i - 1
                  nextj = j + 1
                Else If (dir(i,j)==32) Then
                  nexti = i - 1
                  nextj = j - 1
                Else If (dir(i,j)==1) Then
                  nexti = i
                  nextj = j + 1
                Else If (dir(i,j)==16) Then
                  nexti = i
                  nextj = j - 1
                Else If (dir(i,j)==8) Then
                  nexti = i + 1
                  nextj = j - 1
                Else If (dir(i,j)==4) Then
                  nexti = i + 1
                  nextj = j
                Else If (dir(i,j)==2) Then
                  nexti = i + 1
                  nextj = j + 1
                End If
                i = nexti
                j = nextj
                If (net(nexti,nextj)>0 .And. net(nexti,nextj)/=nbasin(k)) Then
                  basinpnext(k) = net(nexti, nextj)
                  m = 1
                End If

              End Do

            End If

          End Do
        End Subroutine links

        Subroutine leveldetermine(nbasinup, basinup, basinpnext, level, basinlevel, nsub, outlet, drainage, pbasin, numlev)
          Integer, Parameter :: maxlink = 10
          Integer, Parameter :: maxb = 500
          Integer basinpnext(1000) ! next basin number
          Integer nbasinup(1000) ! number of upbasin
          Integer basinup(1000, maxlink) ! up basin number
          Integer basinlevel(1000) ! basin level
        ! total level
          Integer level
          Integer pbasin(1000)
          Integer nsub
          Integer i, j, k, p, q, temp, r
          Integer bs(maxb)
          Integer bn
          Integer ac(maxb)
          Integer drainage(1000)
          Integer w,bm
          Integer tmp(maxb)
          Integer outlet
          Integer numlev(10)
          bn = 0
          bs = -100
          Do k = 1, nsub
            If (basinlevel(k)==1) Then
              bn = bn + 1
              bs(bn) = k
              ac(bn) = drainage(k)
            End If
          End Do

          Do k = 1, bn
            i = bs(k)
            bs(k) = basinpnext(i)
            ac(k) = drainage(bs(k))
          End Do
          Do i = 1, bn - 1
            Do j = i + 1, bn
              If (bs(j)==bs(i) .And. bs(j)>0) bs(j) = -100
            End Do
          End Do
          bm = 0
          Do k = 1, bn
            If (bs(k)>0) Then
              bm = bm + 1
              tmp(bm) = bs(k)
            End If
          End Do
          bn = bm
          Do k = 1, bn
            bs(k) = tmp(k)
            ac(k) = drainage(bs(k))
          End Do
          Call order(ac, bs, bn)

          Do While (bn>0)
            Do k = 1, bn
              i = bs(k)
              If (nbasinup(i)==1) Then
                j = basinup(i, 1)
                basinlevel(i) = basinlevel(j)
              Else If (nbasinup(i)==2) & ! two up basin
                  Then
                j = basinup(i, 1)
                p = basinup(i, 2)
                If (basinlevel(j)==basinlevel(p)) Then
                  If (basinlevel(j)>0) basinlevel(i) = basinlevel(j) + 1
                  If (basinlevel(j)<0) basinlevel(i) = -9999
                Else
                  j = basinup(i, 1)
                  p = basinup(i, 2)
                  If (basinlevel(j)>0 .And. basinlevel(p)>0) Then
                    basinlevel(i) = max(basinlevel(j), basinlevel(p))
                  Else
                    basinlevel(i) = -9999
                  End If
                End If

              Else If (nbasinup(i)==3) & ! three up basin
                  Then
                j = basinup(i, 1)
                p = basinup(i, 2)
                q = basinup(i, 3)
                If (basinlevel(j)==basinlevel(p) .And. basinlevel(p)==basinlevel(q)) Then

                  basinlevel(i) = basinlevel(j) + 1
                Else If (basinlevel(j)==basinlevel(p)) Then
                  If (basinlevel(j)>basinlevel(q)) basinlevel(i) = basinlevel(j) + 1
                  If (basinlevel(j)<basinlevel(q)) basinlevel(i) = basinlevel(q)
                Else If (basinlevel(j)==basinlevel(q)) Then
                  If (basinlevel(j)>basinlevel(p)) basinlevel(i) = basinlevel(j) + 1
                  If (basinlevel(j)<basinlevel(p)) basinlevel(i) = basinlevel(p)
                Else If (basinlevel(p)==basinlevel(q)) Then
                  If (basinlevel(p)>basinlevel(j)) basinlevel(i) = basinlevel(p) + 1
                  If (basinlevel(p)<basinlevel(j)) basinlevel(i) = basinlevel(j)
                Else
                  basinlevel(i) = max(basinlevel(j), basinlevel(p), basinlevel(q))

        !
                End If
              Else If (nbasinup(i)==4) & ! four up basin
                  Then
                j = basinlevel(basinup(i,1))
                p = basinlevel(basinup(i,2))
                q = basinlevel(basinup(i,3))
                r = basinlevel(basinup(i,4))
                If ((j-p)*(p-q)*(q-r)*(j-q)*(j-r)*(p-r)/=0) Then
                  basinlevel(i) = max(j, p, q, r)
                Else If (j==p .And. p==q) Then
                  basinlevel(i) = max(j+1, r)
                Else If (j==q .And. q==r) Then
                  basinlevel(i) = max(j+1, p)
                Else If (p==q .And. q==r) Then
                  basinlevel(i) = max(p+1, j)
                Else If (j==p .And. p==r) Then
                  basinlevel(i) = max(j+1, q)
                Else If (j==p) Then
                  If (q==r) Then
                    basinlevel(i) = max(j+1, q+1)
                  Else
                    basinlevel(i) = max(j+1, q, r)
                  End If
                Else If (j==q) Then
                  If (p==r) Then
                    basinlevel(i) = max(j+1, p+1)
                  Else
                    basinlevel(i) = max(j+1, p, r)
                  End If
                Else If (j==r) Then
                  If (p==q) Then
                    basinlevel(i) = max(j+1, p+1)
                  Else
                    basinlevel(i) = max(j+1, p, q)
                  End If
                Else If (p==q) Then
                  If (j==r) Then
                    basinlevel(i) = max(j+1, p+1)
                  Else
                    basinlevel(i) = max(p+1, j, r)
                  End If
                Else If (p==r) Then
                  If (q==j) Then
                    basinlevel(i) = max(j+1, p+1)
                  Else
                    basinlevel(i) = max(p+1, j, q)
                  End If
                Else If (q==r) Then
                  If (p==j) Then
                    basinlevel(i) = max(j+1, q+1)
                  Else
                    basinlevel(i) = max(q+1, j, p)
                  End If
                End If


              End If
        ! all the upbasin loop
            End Do
            Print *, bn
            w = bn
            If (w==1) Then
              i = bs(bn)
              If (basinpnext(i)<0) Then
                bn = 0
              Else
                bs(bn) = basinpnext(i)
                ac(bn) = drainage(basinpnext(i))
              End If

            Else If (w>1) Then
              Do k = 1, bn
                i = bs(k)
                If (basinpnext(i)>0) Then
                  bs(k) = basinpnext(i)
                  ac(k) = drainage(basinpnext(i))
                Else
                  bs(k) = -100
                  ac(k) = -9999
                End If
              End Do
              bm = 0
              Do i = 1, bn - 1
                Do j = i + 1, bn
                  If (bs(j)==bs(i)) bs(j) = -100

                End Do
              End Do
              Do k = 1, bn
                If (bs(k)>0) Then
                  bm = bm + 1
                  tmp(bm) = bs(k)
                End If
              End Do
              bn = bm
              Do k = 1, bn
                bs(k) = tmp(k)
                ac(k) = drainage(bs(k))
              End Do
              If (bn>1) Call order(ac, bs, bn)
            End If

          End Do
        ! while
          level = basinlevel(outlet)
        !ccc  determine  pbasin number
          Do k = 1, 1
            bn = 0
            Do i = 1, nsub
              If (basinlevel(i)==k) Then
                bn = bn + 1
                pbasin(i) = 1000*k + bn
              End If
            End Do
            Print *, 'basin number of level', k, '  =', bn
            numlev(k) = bn
          End Do
          Do k = 2, level
            bn = 0
            Do i = 1, nsub
              If (basinlevel(i)==k) Then
                bn = bn + 1
                bs(bn) = i
                ac(bn) = drainage(i)
                pbasin(i) = 1000*k + bn
              End If
            End Do
            Print *, 'basin number of level', k, '  =', bn
            numlev(k) = bn
            Do j = 1, bn - 1
              Do p = j + 1, bn
                If (p/=j) Then
                  If (ac(j)>ac(p)) Then
                    q = bs(j)
                    r = bs(p)
                    temp = pbasin(q)
                    pbasin(q) = pbasin(r)
                    pbasin(r) = temp
                    temp = q
                    bs(j) = r
                    bs(p) = temp
                    temp = ac(j)
                    ac(j) = ac(p)
                    ac(p) = temp
                  End If
                End If
        ! p
              End Do
        ! j
            End Do
        ! do 2,level
          End Do
        End Subroutine leveldetermine

        subroutine GetPara()
          implicit none
          Character *5 ncols, nrows
          Character *12 xllcorner, yllcorner
          Character *8 cellsize
          Character *12 nodata
          Real x0, y0
          Integer s, znodata
          Integer, Parameter :: maxlink = 10
          Integer acc(1000, 1000), dir(1000, 1000), net(1000, 1000) !
          Integer subbasin(1000, 1000) ! sub-basin file
          Integer psubbasin(1000, 1000)
          Integer maxtmp
        ! number of sub-basin
          Integer nsub
          Integer nbasin(1000) ! initial number of subbasin
          Integer netbasinstart(1000, 2) ! network start point of each subba
          Integer netbasinend(1000, 2) ! network end point of each subbasi
          Integer i, j, k, tmpi, tmpj,p
          Integer basinpnext(1000) ! next basin number
          Integer drainage(1000) ! drainage basin number
          Integer nbasinup(1000) ! number of upbasin
          Integer basinup(1000, maxlink) ! up basin number
          Integer pbasinup(1000, maxlink) ! up basin pnumber
          Integer pbasin(1000) ! basin pnumber
          Integer basinlevel(1000) ! basin level
        ! total level
          Integer level
          Integer numlev(10)
          Character *12 nstr
          Common /acc/dir
          Integer nc, nr, outlet


          write(*,"('Gis Directory:', A100)")gis_dir
          write(*,"('Parameters Directory:', A100)")para_dir

          Open (1, File=trim(gis_dir)//'acc1k.asc', Status='old') !input1  mqh
          Read (1, *) ncols, nc
          Read (1, *) nrows, nr
          Read (1, *) xllcorner, x0
          Read (1, *) yllcorner, y0
          Read (1, *) cellsize, s
          Read (1, *) nodata, znodata
          Print *, 'x0=,y0=', x0, y0
          Print *, 'nodata,znodata=', nodata, znodata
          Write (*, *)
          Do i = 1, nr
            Read (1, *)(acc(i,j), j=1, nc)
          End Do
          Close (1)

          Open (1, File=trim(gis_dir)//'net1k.asc', Status='old') !input2  mq
          Read (1, *) ncols, nc
          Read (1, *) nrows, nr
          Read (1, *) xllcorner, x0
          Read (1, *) yllcorner, y0
          Read (1, *) cellsize, s
          Read (1, *) nodata, znodata
          Print *, 'x0=,y0=', x0, y0
          Print *, 'nodata,znodata=', nodata, znodata
          Write (*, *)
          Do i = 1, nr
            Read (1, *)(net(i,j), j=1, nc)
          End Do
          Close (1)

          Open (1, File=trim(gis_dir)//'wshed1k.asc', Status='old') !input3
          Read (1, *) ncols, nc
          Read (1, *) nrows, nr
          Read (1, *) xllcorner, x0
          Read (1, *) yllcorner, y0
          Read (1, *) cellsize, s
          Read (1, *) nodata, znodata
          Print *, 'x0=,y0=', x0, y0
          Print *, 'nodata,znodata=', nodata, znodata
          Write (*, *) nc, nr, x0, y0, s
          Do i = 1, nr
            Read (1, *)(subbasin(i,j), j=1, nc)
          End Do
          Close (1)

          Open (1, File=trim(gis_dir)//'dir1k.asc', Status='old') !input4  mq
          Read (1, *) ncols, nc
          Read (1, *) nrows, nr
          Read (1, *) xllcorner, x0
          Read (1, *) yllcorner, y0
          Read (1, *) cellsize, s
          Read (1, *) nodata, znodata
          Print *, 'x0=,y0=', x0, y0
          Print *, 'nodata,znodata=', nodata, znodata
          Write (*, *)
          Do i = 1, nr
            Read (1, *)(dir(i,j), j=1, nc)
          End Do
          Close (1)

        !ccc  rename the net file with subbasin number
          Do i = 1, nr
            Do j = 1, nc
              If (net(i,j)>0) net(i, j) = subbasin(i, j)
            End Do
          End Do

          Open (1, File=trim(gis_dir)//'net1k.asc', Status='replace')
          Write (1, '(A, i16)') ncols, nc
          Write (1, '(A, i16)') nrows, nr
          Write (1, '(A, f24.8)') xllcorner, x0
          Write (1, '(A, f24.8)') yllcorner, y0
          Write (1, '(A, i16)') cellsize, s
          Write (1, '(A, i16)') nodata, znodata
          Do i = 1, nr
            Write (1, *)(net(i,j), j=1, nc)
          End Do
          Close (1)

        !ccc  calculate the number of subbasin
          tmpi = 0
          nsub = 0
          Do i = 1, nr
            Do j = 1, nc
              If (subbasin(i,j)>tmpi) Then
                nsub = subbasin(i, j)
                tmpi = subbasin(i, j)
              End If
            End Do
          End Do

          Print *, 'number of subbasin', nsub
          Do k = 1, nsub
            nbasin(k) = k
          End Do

        !ccc  determine the start and end point of network for each subbasin
          Do k = 1, nsub
            tmpi = nr*nc
            tmpj = -1
            Do i = 1, nr
              Do j = 1, nc
                If (net(i,j)==nbasin(k)) Then
                  If (acc(i,j)>tmpj) Then
                    netbasinend(k, 1) = i
                    netbasinend(k, 2) = j
                    tmpj = acc(i, j)
                  End If
                  If (acc(i,j)<tmpi) Then
                    netbasinstart(k, 1) = i
                    netbasinstart(k, 2) = j
                    tmpi = acc(i, j)
                  End If
        ! first if
                End If
              End Do
            End Do

          End Do
        ! k
          Open (17, File=trim(para_dir)//'basin_network.txt')
          Do k = 1, nsub
            Write (17, *) k, netbasinend(k, 1), netbasinend(k, 2)
          End Do
          Close (17)
          maxtmp = 0
          Do k = 1, nsub
            i = netbasinend(k, 1)
            j = netbasinend(k, 2)
            drainage(k) = acc(i, j)
            If (acc(i,j)>maxtmp) Then
              maxtmp = acc(i, j)
              outlet = k
            End If
          End Do
          Print *, drainage(outlet), outlet

        !ccc  detemine the link of each subbasin
          Call links(dir, netbasinend,net, nsub, nbasin, basinpnext, outlet)
          nbasinup = 0
          basinup = -9999
          basinlevel = -9999
          pbasinup = -9999
          Do k = 1, nsub
            Do i = 1, nsub
              If (basinpnext(i)==k) Then
                nbasinup(k) = nbasinup(k) + 1
                basinup(k, nbasinup(k)) = i
              End If
            End Do

          End Do
        ! k
          Do k = 1, nsub, 1
            If (nbasinup(k)==0) Then
              basinlevel(k) = 1
            End If

          End Do
        ! k
          Open (18, File=trim(para_dir)//'basinup.txt', Status='replace')
          Do k = 1, nsub
            Write (18, *) nbasin(k), nbasinup(k)
            If (nbasinup(k)>0) Then
              Do i = 1, nbasinup(k)
                Write (18, *) basinup(k, i)
              End Do
            End If
          End Do
          Close (18)

        !ccc  determine the level
          Call leveldetermine(nbasinup, basinup, basinpnext, level, basinlevel,nsub, outlet, drainage, pbasin, numlev)
          pbasinup = -9999
          Do k = 1, nsub
            If (nbasinup(k)>0) Then
              Do i = 1, nbasinup(k)
                j = basinup(k, i)
                pbasinup(k, i) = pbasin(j)
              End Do
            End If
          End Do

          Open (18, File=trim(para_dir)//'basinlevel.txt', Status='replace')
          Do k = 1, nsub
            Write (18, *) nbasin(k), basinlevel(k), pbasin(k)
          End Do
          Close (18)

          psubbasin = -9999
          Do i = 1, nr
            Do j = 1, nc
              Do k = 1, nsub
                If (subbasin(i,j)==nbasin(k)) Then
                  psubbasin(i, j) = pbasin(k)
                End If
              End Do
            End Do
          End Do

        !ccc  out put
          Open (1, File=trim(para_dir)//'pbasin.asc', Status='replace')
          Write (1, '(A, i16)') ncols, nc
          Write (1, '(A, i16)') nrows, nr
          Write (1, '(A, f24.8)') xllcorner, x0
          Write (1, '(A, f24.8)') yllcorner, y0
          Write (1, '(A, i16)') cellsize, s
          Write (1, '(A, i16)') nodata, znodata

          Do i = 1, nr
            Write (1, *)(psubbasin(i,j), j=1, nc)
          End Do
          Close (1)

          j = 0
          Open (1, File=trim(para_dir)//'subcatchment.txt', Status='replace')
          Do k = 1, level
            Do i = 1, numlev(k)
              j = j + 1
              nstr = 'ws'
              Write (1, '(i4,3x,a2,i4,3x,i4.4)') j, nstr, k*1000 + i, k*1000 + i
            End Do
          End Do
          Close (1)
          j = 0
          Open (1, File=trim(para_dir)//'subbasin.dat', Status='replace')
          Do k = 1, level
            Do i = 1, numlev(k)
              j = j + 1
              Do p = 1, nsub
                If (pbasin(p)==k*1000+i) Then
                  If (nbasinup(p)==0) Then
                    Write (1, '(i4,2x,i4,2x,4i4)') j, k*1000 + i, 0, 0, 0, 0
                  Else If (nbasinup(p)==1) Then
                    Write (1, '(i4,2x,i4,2x,4i4)') j, k*1000 + i, pbasinup(p, 1), 0, 0, 0
                  Else If (nbasinup(p)==2) Then
                    Write (1, '(i4,2x,i4,2x,4i6)') j, k*1000 + i, pbasinup(p, 1), pbasinup(p, 2), 0, 0
                  Else If (nbasinup(p)==3) Then
                    Write (1, '(i4,2x,i4,2x,4i6)') j, k*1000 + i, pbasinup(p, 1), pbasinup(p, 2), pbasinup(p, 3), 0
                  Else If (nbasinup(p)==4) Then
                    Write (1, '(i4,2x,i4,2x,4i6)') j, k*1000 + i, pbasinup(p, 1), pbasinup(p, 2), pbasinup(p, 3), pbasinup(p, 4)
                  End If
                End If
              End Do
            End Do
          End Do
          Close (1)

          Open (1, File=trim(para_dir)//'drainage.txt', Status='replace')
          Do k = 1, level
            Do i = 1, numlev(k)
              j = j + 1
              Do p = 1, nsub
                If (pbasin(p)==k*1000+i) Then
                  Write (1, *) pbasin(p), drainage(p)
                End If
              End Do
            End Do
          End Do
          Close (1)
        End subroutine
end module Priver


subroutine morph(nsub,dem_home,model_dir)
  implicit none
  integer,parameter::ntl=5
  integer nsub,num_sub
  integer ia1,ia2

  Character *200 dem_home, model_dir
  Character *6 ws_name(nsub)

  Character *200 rivername
  Character *200 dir1, dir2
  Integer i, j

  Common /river/rivername, dir1, dir2

  Integer kfs1(10), kfs2(10, 10), kfs3(10, 10, 10)
  Data kfs1/10*0/
  Data kfs2/100*0/
  Data kfs3/1000*0/

!                 *** CHANGE THIS PART ***
!-----------------------------------------------------------------------
!  home = '..\' !KEEP in LOCAL dir
!  dem_home = '..\subs\' !KEEP in LOCAL dir
!  model_dir = '..\' !-----------------------------------------------------------------------
!     generate subcatchment name from kfs.dat
!KEEP in LOCAL dir
  Call strlen(model_dir, ia1, ia2)
  Open (1, File=model_dir(ia1:ia2)//'subcatchment.txt', Status='old')
  Do num_sub = 1, nsub
    Read (1, *) i, ws_name(num_sub), j
  End Do
  Close (1)
!
!-----------------------------------------------------------------------
!
  dir2 = model_dir
  Do i = 1, nsub
    rivername = ws_name(i)
    Print *, i, '     ', rivername(1:6)
    Call strlen(dem_home, ia1, ia2)
    dir1 = dem_home(ia1:ia2) // 'sub_' // ws_name(i) // '/'
    Call basin_morph
  End Do
End subroutine morph
!
!************************************************************************
Subroutine basin_morph
  Character *5 ncols, nrows
  Character *9 xllcorner, yllcorner
  Character *8 cellsize
  Character *12 nodata
  Real distmp(800, 800), bedslope_tmp(800, 800)
  Integer location_tmp(800, 800), location(500000)
  Integer row(500000), col(500000)
  Real dis(500000), bedslope(500000)

  Character *200 infile(3)
  Character *200 outfile, rivername
  Character *200 dir1, dir2

  Common /river/rivername, dir1, dir2

  infile(1) = 'distance.asc'
  infile(2) = 'watershed.asc'
  infile(3) = 'bedslope.asc'
  Call strlen(rivername, ia1, ia2)
  outfile = rivername(ia1:ia2) // '_morph'

  Call strlen(dir2, ib1, ib2)
  Open (3, File=dir2(ib1:ib2)//outfile, Status='unknown')
  Do k = 1, 3
    Call strlen(dir1, ib1, ib2)
    Open (1, File=dir1(ib1:ib2)//infile(k), Status='old')
    Read (1, '(A, i16)') ncols, nc
    Read (1, '(A, i16)') nrows, nr
    Read (1, '(A, f16.3)') xllcorner, x0
    Read (1, '(A, f16.3)') yllcorner, y0
    Read (1, '(A, f16.0)') cellsize, s
    Read (1, '(A, f16.0)') nodata, znodata
!
!     !CHECK NCOLS AND NROWS HERE!

    Write (*, *) 'nr=,nc=', nr, nc

    Do i = 1, nr
      If (k==1) Read (1, *)(distmp(i,j), j=1, nc)
!      write(*,*)  (distmp(i,j), j=nc,nc)
      If (k==2) Read (1, *)(location_tmp(i,j), j=1, nc)
!      write(*,*)  (location_tmp(i,j), j=nc,nc)
      If (k==3) Read (1, *)(bedslope_tmp(i,j), j=1, nc)
!      write(*,*)  (bedslope_tmp(i,j), j=nc,nc)
!      pause
    End Do
    Close (1)
  End Do
  m = 0

  Write (*, *) 'nr=,nc=', nr, nc
!      pause

  Do i = 1, nr
    Do j = 1, nc

      If (location_tmp(i,j)/=-9999 .And. (distmp(i,j)==-9999 .Or. bedslope_tmp(i,j)==-9999)) Then
        Print *, 'wrong in basin_morph:', rivername, i, j, distmp(i, j), location_tmp(i, j), bedslope_tmp(i, j)
        Stop
      End If

      If (location_tmp(i,j)/=-9999) Then
        m = m + 1
        dis(m) = distmp(i, j)
        location(m) = location_tmp(i, j)
        row(m) = i
        col(m) = j
        bedslope(m) = bedslope_tmp(i, j)

        Write (3, 10) dis(m), row(m), col(m), bedslope(m)
      End If
    End Do
  End Do
  Close (3)
  Return
  10 Format (F15.3, 2I8, F15.5)
End Subroutine basin_morph
!
!     sub-program to detect the string length
Subroutine strlen(str, l1, l2)
  Character *200 str
  Integer i, l1, l2, k
  k = 0
  Do i = 1, 200
    If (k==0 .And. str(i:i)/=' ') Then
      l1 = i
      k = 1
    Else If (k==1 .And. str(i:i)==' ') Then
      l2 = i - 1
      Return
    End If
  End Do
  l2 = i
  Return
End Subroutine strlen


module parameters
    implicit none
    integer::nc=100,nr=100,nsub=50
    character*200::home="./"
    Character*100 subcatch
    real dx_max,width_min, width_max, height_min, height_max, roughness_min, roughness_max

    contains

        subroutine para()
          implicit none
          integer ia1,ia2,ii
          integer,parameter::ntl=5
          Character *6 ws_name(nsub)
        !
        !     ** for 1-km grid, RIVER LENGTH range from 5x to 10x **
          dx_max = 5000.0
        !----------------------------------------------------------------------
        !     specify river cross-section parameters
              !para_dir           =      "../parameter/"
        ! home = '../' !!KEEP in LOCAL dir!! (change3£©
          Call strlen(home, ia1, ia2)
          Print *, home(ia1:ia2)
          Open (110, File=home(ia1:ia2)//'subcatchment.dat', Status='old')

          Read (110, *)
          10 Read (110, *, End=100) ii, ws_name(ii), width_min, width_max, height_min, height_max, roughness_min, roughness_max

          subcatch = ws_name(ii)
          Call parameter_model(ii)
          Goto 10
          100 Close (110)
        End subroutine
        !
        !********************************************************************
        !
        Subroutine parameter_model(ii)
          implicit none
          integer ii

          Print *
          Print *, ii, '   ', subcatch(1:10)
          Print *, width_min, width_max, height_min, height_max, roughness_min, roughness_max
        !
          Print *, 'Calculating sub-basin parameters...'
          Print *
          Call basin
          Return
        End Subroutine parameter_model
        !**********************************************************************
        !                                                                     *
        !                 ***** Basin Parameter *****                         *
        !                                                                     *
        !**********************************************************************
        !   Sub-program for calculating the sub-basin parameters using the
        !   geomorphology data ('##_morph') generated by GIS.
        !
        !   Variables:
        !   dx, s0:          length and average slope for each flow
        !                     interval in simulation order(outlet-->source)

        !   ddx, ss0:        length, average slope of local river interval
        !                     in area function order(source-->outlet)

        !   num:             number of flow intervals for hydrologic model
        !   distance:        total distance of the river (equal to maxmum river
        !                   length: dis_max).
        !
        !**********************************************************************
        !
        Subroutine basin
          implicit none
          Character *100 output, basinmorph
          Real dis(500000), bed_angle(500000)
          Integer ir(500000), ic(500000)
          Real dx(1000), s0(1000), b(1000), dr(1000), roughness(1000)
          Integer np(1000), ipr(1000, 50000), ipc(1000, 50000)

        !
          real d1,d2,dis_max,distance,dtmp,tmp_angle
          integer i,ia1,ia2,iarea,ib1,ib2,j,ii,n_totalgrid,num

              !NEED CHANGE HERE!
        ! dir = '../' !KEEP in LOCAL dir
          Call strlen(home, ia1, ia2)

          Call strlen(subcatch, ib1, ib2)
          basinmorph = home(ia1:ia2) // subcatch(ib1:ib2) // '_morph'
          Print *, trim(subcatch) // '_morph'

          output = home(ia1:ia2) // subcatch(ib1:ib2) // '_river'
        !
          dis_max = 0.0
          Open (1, File=basinmorph, Status='old')
          i = 1
          10 Read (1, *, End=20) dis(i), ir(i), ic(i), bed_angle(i)
        !
              !NEED CHANGE HERE!
          If (ir(i)<=0 .Or. ir(i)>nr .Or. ic(i)<=0 .Or. ic(i)>nc) Then
            Print *, 'wrong grid number:', i, ir(i), ic(i)
            print *, 'Change nc and nr setting.'
            Stop
          End If
          If (dis(i)>dis_max) dis_max = dis(i)
          i = i + 1
          Goto 10
          20 Close (1)
          n_totalgrid = i - 1

          Print *, 'n_totalgrid=', n_totalgrid
        !
        !     calculate river interval parameters (interval-length <= dx_max)
          num = int(dis_max/dx_max) + 1

          d1 = 0.0
          d2 = 0.0
          Do i = 1, num
            ii = num - i + 1
            dx(ii) = dx_max
            d2 = d1 + dx(ii)
            If (i==num) d2 = d1 + dx(ii) + 1.0
            If (i==num) dx(ii) = dis_max - d1
            np(ii) = 0
            tmp_angle = 0.0
            Do j = 1, n_totalgrid
              If (dis(j)>=d1 .And. dis(j)<d2) Then
                np(ii) = np(ii) + 1
                ipr(ii, np(ii)) = ir(j)
                ipc(ii, np(ii)) = ic(j)
                If (ir(j)<=0 .Or. ic(j)<=0) Print *, 'wrong:', ii, j, ir(j), ic(j)
                tmp_angle = tmp_angle + bed_angle(j)
              End If
            End Do
            iarea = iarea + np(ii)
            If (np(ii)==0) Then
              Print *, ii, np(ii)
              np(ii) = 1
            End If
            If (dx(ii)==0) Then

              dx(ii) = 500
            End If
            tmp_angle = tmp_angle/float(np(ii))
            s0(ii) = tmp_angle/100.0
            d1 = d2
          End Do
        !
          distance = 0.0
          Do i = 1, num
            distance = distance + dx(i)
            If (s0(i)<0.0) Print *, 'Error in calculating river bed slope', i, s0(i)
          End Do
        !
          If (abs(distance-dis_max)>1.0) Print *, 'Error in calculating main channel distance IN RIVER', distance, dis_max

          dtmp = 0.0
          Do i = 1, num
            dtmp = dtmp + dx(i)
            b(i) = width_max + float(num-i)*(width_min-width_max)/float(num)
            roughness(i) = roughness_min + float(num-i)*(roughness_max-roughness_min)/float(num)
            dr(i) = height_max + float(num-i)*(height_min-height_max)/float(num)
        !

            If (b(i)<=0. .Or. roughness(i)<=0. .Or. dr(i)<=0.0) Print *, 'wrong:  ', subcatch(ia1:ia2)
          End Do
        !
          Call strlen(output, ia1, ia2)
          Open (3, File=output(ia1:ia2), Status='unknown')
          Write (3, *) num
          Do i = 1, num
            Write (3, 200) np(i), dx(i), s0(i), b(i), roughness(i), dr(i)
            Write (3, *)(ipr(i,j), ipc(i,j), j=1, np(i))
          End Do
          Close (3)
          Return
          200 Format (I5, F12.2, F12.6, F10.2, F10.5, F10.2)
        End Subroutine basin


end module

program main
    use parameters
    Character *200 dem_home, model_dir
    dem_home="D:\pythowork\workspace1\station\station0\morphx\subs"
    model_dir="D:\pythowork\workspace1\station\station0\morphx\"
    call morph(17,dem_home,model_dir)
end program
