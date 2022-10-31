        subroutine gaussfilter1(mat,mx,my,mz,mw,filtermat,ksiz,mask,smoo
     1  thed) 

COMMENT applies a 3D discrete filter to a 4D array using a mask

        integer mx,my,mz,mw,ksiz,m1,m2,m3,e,i,j,k,m
        double precision mat(mx,my,mz,mw),filtermat(ksiz,ksiz,ksiz)
        double precision smoothed(mx,my,mz,mw),mask(mx,my,mz),total,f


        e=(ksiz+1)/2
        do i=1,mx
           do j=1,my
              do k=1,mz
                 if(mask(i,j,k).eq.1.0) then
                    do m=1,mw
                       total=0.0
                       f=0.0

                       do m1=1,ksiz
                          do m2=1,ksiz
                             do m3=1,ksiz

       if((0.lt.(i-e+m1)).and.((i-e+m1).le.mx).and.(0.lt.(j-e+m2))) then
       if(((j-e+m2).le.my).and.(0.lt.(k-e+m3)).and.((k-e+m3).le.mz))then


        total=total+filtermat(m1,m2,m3)*mat(i-e+m1,j-e+m2,k-e+m3,m)*
     1  mask(i-e+m1,j-e+m2,k-e+m3)

        f=f+filtermat(m1,m2,m3)*mask(i-e+m1,j-e+m2,k-e+m3)

                                   end if
                                end if
                             end do
                          end do
                       end do

                       smoothed(i,j,k,m)=(total*mask(i,j,k))/f

                    end do
                 else
                    do m=1,mw
                       smoothed(i,j,k,m)=0.0 
                    end do
                 end if
              end do
           end do
        end do

        return
        end




        subroutine gaussfilter2(mat, mx, my, mz, mw, filtermat, d, 
     1  mask, smoothed)

COMMENT applies a 3D discrete filter to a 4D array and normalises so variance is unchanged
        integer mx, my, mz, mw, d, m1, m2, m3, e, i, j, k, l
        double precision mat(mx, my, mz, mw), filtermat(d, d, d)
        double precision smoothed(mx, my, mz, mw), mask(mx, my, mz),
     1  total, f, f2

        e = (d + 1) / 2
        do i = 1, mx
           do j = 1, my
              do k = 1, mz
                 if(mask(i, j, k).eq.1.0) then
                    do l = 1, mw

                       total = 0.0
                       f = 0.0
                       f2 = 0.0

                       do m1 = 1, d
                          do m2 = 1, d
                             do m3 = 1, d

        if((0.lt.(i - e + m1)).and.((i - e + m1).le.mx).and.
     1  (0.lt.(j - e + m2))) then
        if(((j - e + m2).le.my).and.(0.lt.(k - e + m3)).and.
     1  ((k - e + m3).le.mz))then
        total = total + filtermat(m1, m2, m3) * 
     1  mat(i - e + m1, j - e + m2, k - e + m3, l)
COMMENT              f = f + filtermat(m1, m2, m3)
         f2 = f2 + filtermat(m1, m2, m3) ** (2.0)

                                   end if
                                end if
                             end do
                          end do
                       end do

                       smoothed(i, j, k, l) = total / (sqrt(f2))
                    end do
                 else 
                    do l = 1, mw
                       smoothed(i, j, k, l) = 0.0 
                    end do
                 end if
              end do 
           end do
        end do
        
        return
        
        end
