program sc_fermi
!
! Program to determine the self consistent Fermi level in a 
! semiconductor (or metal) by considering the formation energy
! of a set of defects in different charge states and the DOS
! of the perfect system.
!
! The relevant equations are Eq. (2-4) of Yang et al., Scientific
! Reports, 5:16977 (2015).
!
! ************************
! UPDATE - equations in that reference are not correct for n0 and p0
! ************************
!
! The formation energies of the defects in the different charge 
! states are read in at E_F = 0, and the transition levels are 
! determined to calculate the lowest formation energy charge state
! at each value of E_F. 
!
! The concentration of defect alpha in charge q is determined by:
!
! n(alpha,q) = N0 g e^(-H_f(E_F)/kT)
!
! where N is the density of sites where alpha can form and g is the
! degeneracy of the defect, which includes electronic and structural
! degeneracy.
!
! ************************
! NEW ADDITION: use Fermi Dirac statistics - first calculate the
! concentration of the neutral defect:
!
! c_X^0 = N_X g_0 e^(-E_f^0/kT)
!
! Then get the concentrations of the charged defects as:
! (see labbook 24/5/17)
!
! c_X^q = N_X g_q e^(-E_f^q/kT)
!
! N_X is the density of sites where X can form, c_X^q, g_q, and 
! E_f^q are the concentration, degeneracy factor and formation 
! energy of X in charge state q
!
! (If E_f^0 is not supplied, set C_X^0 = 0)
!
! This analysis follows from (Ashcroft and Mermin Chapter 28) setting
! the number of electrons introduced by defects using Gibb's factors:
!
! <n> = sum_j N_j exp-(E_j - E_F N_j)/kT / sum_j exp-(E_j - E_F N_j)/kT
!
! with N_j = q, E_j = -qe(q|0)
!
! The concentration of defect X is then:
!
! C_X = C_X^0 + sum_q C_X^q
!
! Return to the user the total, and the concentration per charge state
! ************************
!
! Electron and hole concentrations n0 and p0 are determined by:
!
! n0 = N_c e^(-(Eg-E_F)/kT); N_c = int_Eg^inf dE D(E) / (1+(e^(e-Eg)/kT)
!
! p0 = N_v e^(-(Eg-E_F)/kT); N_v = int_-inf^0 dE D(E) / (1+(e^(-E)/kT)
!
! ************************
! ABOVE INCORRECT - the equations are as 28.9 in Ashcroft Mermin
!
! n0 = int_Eg^inf dE D(E) / (1+e^(E-E_F)/kT)
!
! p0 = int_-inf^0 dE D(E) / (1+e^(E_F-E)/kT)
!
! where D(E) is the density of states of the perfect system.
! ************************
!
! The self-consistent E_F is determined where 
!
! n0 + n(alpha,q)[acceptors]|q| = p0 + n(alpha,q)[donors]|q|  (*)
!
! by a scanning technique.
!
! The H_f(alpha,q) at E_F=0, N0, D(E), g, and T are read in, while
! the SC E_F and concentrations are output.
! 
! ************************
! UPDATE Aug 30th 2019: as pointed out by Ben Morgan, there was a 
! bug in the DOS integrations, where the summations up to the Fermi
! level cut off one slice short of the correct amount. This bug
! has now been fixed.
!
! J. Buckeridge June 2016
!
! UPDATE - May 2017
! UPDATE - Apr 2018
! BUGFIX - Aug 2019
!
  implicit none
!
  integer, parameter :: nmax = 25000
  integer, parameter :: nless = 500
  integer, parameter :: bigint = 1e6
  real*8, parameter :: small = 1.d-14
  real*8, parameter :: lsmall = 1.d-7
  real*8, parameter :: vlarge = 1.014232054735d304
  real*8, parameter :: largefac = 700.d0
  real*8, parameter :: kboltz = 8.617343d-5 ! Boltzmann constant (eV/K)
  integer i, j, k, l, stat, stat1, stat2
  integer count, lownum, kcount, transnum(nless), linelength, pspace
  integer ndefects, nspinpol, ncharge(nless), nsite(nless), numdos
  integer kmax, kmin
  integer countdos
  real*8 nelect, egap, charge(nless,nless), energy(nless,nless), deg(nless,nless)
  real*8 temperature, ener0(nless), deg0(nless), enerq
  real*8 pi, sum1, temp1, emin, emax, loweng, fermitrans(nless), enertmp
  real*8 lattvec(3,3), factor, volume
  real*8 edos(nmax), dos(nmax), dosu(nmax), dosd(nmax), de, expfac
  real*8 delphi, estep, direction, n0, p0, lhs, rhs, diff, diffold
  real*8 delphi_min, delphibig, dirbig
  real*8 conc0(nless), concd(nless,nless), concall(nless)
  character*3 str1
  character*10 name(nless)
  logical check0, check_under, check_under2
  logical checkconv
  logical checkposcar
  character*80 line
!
  pi = acos(-1.0)
!
! Write header
!
  write(*,'(a)') "**************************************************************"
  write(*,'(a)') 
  write(*,'(a)') "   SSSS    CCCC      FFFFFF  EEEEEE   RRRR   MM     MM  IIIII" 
  write(*,'(a)') "  SS   S  CC   C     FF      EE      RR   R  MMM   MMM    I" 
  write(*,'(a)') "  SS      CC         FF      EE      RR  R   M MM MM M    I" 
  write(*,'(a)') "   SSSS   CC     --- FFFFFF  EEEEEE  RRRR    M  MMM  M    I" 
  write(*,'(a)') "      SS  CC         FF      EE      R   R   M   M   M    I" 
  write(*,'(a)') "  S   SS  CC   C     FF      EE      R   RR  M   M   M    I" 
  write(*,'(a)') "   SSSS    CCCC      FF      EEEEEE  R   RR  M   M   M  IIIII" 
  write(*,'(a)') 
  write(*,'(a)') "Energies in eV, temperature in Kelvin, DOS in states/unit cell" 
  write(*,'(a)') 
  write(*,'(a)') "------"
  write(*,'(a)') "j.buckeridge@ucl.ac.uk 2018"
  write(*,'(a)') "**************************************************************"
  write(*,'(a)') 
!
! Open unitcell.dat to read in vectors and determine volume of pure cell
!
  checkposcar = .false.
  open(unit=11,file="unitcell.dat",status="old",iostat=stat)
!
  if(stat /= 0) then
     write(*,'(a)') "ERROR: unitcell.dat (unit cell for which DOS determined) is needed!!"
     checkposcar = .true.
  else
     write(*,'(a)') "unitcell.dat found..."
     write(*,'(a)') "(Should be cell for which DOS was determined!)"
     write(*,*)
!
! Read in scaling factor from the file
!
     count = 0
     do
        count = count + 1
        read(11,*) line
        if(line(1:1) == "#") cycle
        read(line,*,iostat=stat) factor
        if(stat /= 0) then
           write(*,'(a)') "ERROR: scaling factor for lattice vectors not found!!"
           stop
        else
           exit
        endif
        if(count > bigint) then
           write(*,'(a)') "ERROR: cannot read data in unitcell.dat!!"
           stop
        endif
     enddo
!
! Now read in the lattice vectors
!
     count = 0
     i = 0
     do
        count = count + 1
        read(11,*) line
        if(line(1:1) == "#") then
           cycle
        else
           backspace(unit=11)
           i = i+1
           read(11,*,iostat=stat) (lattvec(i,j),j=1,3)
           if(stat /= 0) then
              write(*,'(a39,x,i1,x,a2)') "ERROR: could not read in lattice vector", i, "!!"
              stop
           endif
           if(i == 3) exit
        endif
        if(count > bigint) then
           write(*,'(a)') "ERROR: lattice vectors not found in unitcell.dat!!"
           stop
        endif
     enddo
     lattvec = lattvec * factor
     close(unit=11)
  endif

!
! Open POSCAR to read in vectors and determine volume of pure cell
!
!
  if(checkposcar) then
     open(unit=11,file="POSCAR",status="old",iostat=stat)
!
     if(stat /= 0) then
        write(*,'(a)') "Exiting as no unit cell data can be found..."
        stop
     endif
     write(*,'(a)') "...but found a POSCAR! Will use that instead..."
     write(*,'(a)') "(Should be cell for which DOS was determined!)"
     write(*,*)
!
! Read in lattice vectors
!
     read(11,*)
     read(11,*) factor
     do i=1,3
        read(11,*) (lattvec(i,j),j=1,3)
     enddo
     lattvec = lattvec * factor
     close(unit=11)
  endif
!
! Calculate cell volume
!
  volume = lattvec(1,1) * (lattvec(2,2)*lattvec(3,3) - lattvec(2,3)*lattvec(3,2))&
       & + lattvec(1,2) * (lattvec(2,3)*lattvec(3,1) - lattvec(2,1)*lattvec(3,3))&
       & + lattvec(1,3) * (lattvec(2,1)*lattvec(3,2) - lattvec(2,2)*lattvec(3,1))
  write(*,'(a15,2x,f13.6,x,a3)') "Volume of cell:", volume, "A^3"
  write(*,*)
!
! Open input file and read in data
!
  open(unit=11, file="input-fermi.dat", status="old", iostat=stat)
  if(stat /= 0) then
     write(*,'(a)') "No input-fermi.dat file found - bye!!"
     stop
  endif
!
! Read in whether system is spin polarised or not (will be used 
! for reading in DOS)
!
  count = 0
  do
     count = count + 1
     read(11,*) line
     if(line(1:1) == "#") cycle
     read(line,*,iostat=stat) nspinpol
     if(stat /= 0) then
        write(*,'(a)') "ERROR: spin polarisation specification not found!!"
        stop
     else
        if(nspinpol == 1) then
           write(*,'(a)') "Found non-spin polarised system..."
           exit
        elseif(nspinpol == 2) then
           write(*,'(a)') "Found spin polarised system..."
           exit
        else
           write(*,'(a)') "ERROR: spin polarisation specification must by 1 or 2!!"
           stop
        endif
     endif
     if(count > bigint) then
        write(*,'(a)') "ERROR: spin polarisation specification not found!!"
        stop
     endif
  enddo
  write(*,*)
!
! Read in number of electrons (will be used to renormalise DOS)
!
  count = 0
  do
     count = count + 1
     read(11,*) line
     if(line(1:1) == "#") cycle
     read(line,*,iostat=stat) nelect
     if(stat /= 0) then
        write(*,'(a)') "ERROR: number of electrons not found!!"
        stop
     else
        write(*,'(a30,2x,f13.6)') "Number of electrons in system:", nelect
        exit
     endif
     if(count > bigint) then
        write(*,'(a)') "ERROR: number of electrons not found!!"
        stop
     endif
  enddo
!
! Read in energy gap
!
  count = 0
  do
     count = count + 1
     read(11,*) line
     if(line(1:1) == "#") cycle
     read(line,*,iostat=stat) egap
     if(stat /= 0) then
        write(*,'(a)') "ERROR: energy gap not found!!"
        stop
     else
        write(*,'(a30,2x,f13.6,x,a2)') "Energy gap of system:", egap, "eV"
        if(egap < 0.d0) then
           write(*,'(a)') "You have a negative gap - this is going to get weird..."
        endif
        exit
     endif
     if(count > bigint) then
        write(*,'(a)') "ERROR: energy gap not found!!"
        stop
     endif
  enddo
!
! Read in temperature
!
  count = 0
  do
     count = count + 1
     read(11,*) line
     if(line(1:1) == "#") cycle
     read(line,*,iostat=stat) temperature
     if(stat /= 0) then
        write(*,'(a)') "ERROR: temperature not found!!"
        stop
     else
        write(*,'(a30,2x,f13.6,x,a1)') "Temperature :", temperature, "K"
        exit
     endif
     if(count > bigint) then
        write(*,'(a)') "ERROR: temperature not found!!"
        stop
     endif
  enddo
!
! Read in number of defect species
!
  count = 0
  do
     count = count + 1
     read(11,*) line
     if(line(1:1) == "#") cycle
     read(line,*,iostat=stat) ndefects
     if(stat /= 0) then
        write(*,'(a)') "ERROR: number of defect species not found!!"
        stop
     else
        write(*,'(a30,2x,i4)') "Number of defect species:", ndefects
        exit
     endif
     if(count > bigint) then
        write(*,'(a)') "ERROR: number of defect species not found!!"
        stop
     endif
     if(temperature <= small) then
        write(*,'(a)') "ERROR: zero temperature means zero defects!!"
        stop
     endif
  enddo
  write(*,*)
!
! For each species read in name, number of sites in the unit cell, 
! charge states, degeneracy and energies
!
  if(ndefects > 0) then
     count = 0
     i = 1
     do
        count = count + 1
        read(11,'(a)') line
        if(line(1:1) == "#") cycle
        read(line,*,iostat=stat) name(i), ncharge(i), nsite(i)
        if(stat /= 0) then
           write(*,'(a)') "ERROR: failed to read in defect name/charge!!"
           stop
        else
           if(ncharge(i) <= 0) then
              write(*,'(a13,x,i4,x,a36)') "ERROR: defect", i, &
                   &"has idiotic number of charge states!!"
              stop
           endif
           j = 1
           do
              read(11,'(a)') line
              if(line(1:1) == "#") cycle
              read(line,*,iostat=stat) charge(i,j), energy(i,j), deg(i,j)
              if(stat /= 0) then
                 write(*,'(a29,i3,x,a4,i3,a2)') "ERROR: failed read for defect", &
                      &i, "line", j, "!!"
                 stop
              endif
              j = j+1
              if(j > ncharge(i)) exit
           enddo
        endif
        if(i >= ndefects) exit
        i = i+1
        if(count > bigint) then
           write(*,'(a41,i3,a2)') "ERROR: failed to find all data for defect", i, "!!"
           stop
        endif
     enddo
!
! For sanity check return the defects and charge states found
!
     write(*,'(a)') "Defects found:"
     write(*,'(a4,2x,a15,2x,a20)') "Name", "# charge states", "# sites in unit cell"
     do i=1,ndefects
        write(*,'(a10,7x,i4,18x,i4)') name(i), ncharge(i), nsite(i)
     enddo
     write(*,*)
  else
     write(*,'(a)') "0 defects found..."
     write(*,*)
  endif
  close(unit=11)
!
! Read in density of states and renormalise. The energy scale should 
! be relative to the VBM (or Fermi level for a metal). Take into account
! spin polarisation specification
!
  open(unit=11, file="totdos.dat", status="old", iostat=stat)
  if(stat /= 0) then
     write(*,'(a)') "ERROR: no totdos.dat file found - bye!!"
     stop
  endif
!
  i = 0
  countdos = 0
  do
     read(11,'(a)',iostat=stat) line
     if(stat < 0) exit
     if(line(1:1) == "#") cycle
     i = i+1
     if(nspinpol == 1) then
        read(line,*) edos(i), dos(i)
        if(dos(i) < 0.d0) then
           countdos = countdos + 1
        endif
     else
        read(line,*) edos(i), dosu(i), dosd(i)
        if(dosu(i) < 0.d0 .or. dosd(i) < 0.d0) then
           countdos = countdos + 1
        endif
        dos(i) = dosu(i) + dosd(i)
     endif
  enddo
  close(unit=11)
  if(countdos > 0) then
     write(*,'(a)') "WARNING: Found negative value(s) of DOS!"
     write(*,'(a)') "         Do you know what you are doing?"
     write(*,'(a)') "         These may cause serious problems..."
     write(*,*)
  endif
  numdos = i
  emin = edos(1)
  emax = edos(numdos)
  if(egap > emax) then
     write(*,'(a)') "ERROR: your conduction band is not present in energy range of DOS!!"
     stop
  endif
!
! For each defect calculate the transition levels and write them to 
! file for plotting. Determine them by finding the lowest energy defect
! at emin, then find earliest transition to a less positive charge state,
! and continue until all transitions are found
!
  if(ndefects > 0) then
     open(unit=12, file="transition-levels.dat", status="replace")
     do i=1,ndefects
        write(12,'(a1,x,a10)') "#", name(i)
        loweng = energy(i,1) + charge(i,1) * emin
        lownum = 1
        if(ncharge(i) > 1) then
           do j=1,ncharge(i)
              enertmp = energy(i,j) + charge(i,j) * emin
              if(enertmp - loweng <= small) then
                 loweng = enertmp
                 lownum = j
              endif
           enddo
        endif
        write(12,*) emin, loweng
        do
           count = 0
           do j=1,ncharge(i)
              if(j == lownum) cycle
              if(charge(i,j) < charge(i,lownum) .and. &
                   &abs(charge(i,j) - charge(i,lownum)) > small) then
                 count = count + 1
                 fermitrans(count) = (energy(i,j) - energy(i,lownum)) /&
                      &(charge(i,lownum) - charge(i,j))
                 transnum(count) = j
              endif
           enddo
           if(count == 0) exit
           loweng = fermitrans(1)
           kcount = 1
           do k=1,count
              if(fermitrans(k) - loweng <= small) then
                 loweng = fermitrans(k)
                 kcount = k
              endif
           enddo
           write(12,*) fermitrans(kcount), energy(i,lownum) + &
                &charge(i,lownum) * fermitrans(kcount)
           lownum = transnum(kcount)
        enddo
        write(12,*) emax, energy(i,lownum) + charge(i,lownum) * emax
        write(12,*)
     enddo
     close(unit=12)
  endif
!
! Now integrate DOS up to VBM (using trapezoidal rule) to get number of 
! states per unit cell. Renormalise according to number of electrons
!
  sum1 = 0.d0
  i=0
  do
     i = i+1
     if(edos(i+1) > 0.d0) exit
     de = (edos(i+1) - edos(i))
     sum1 = sum1  + (dos(i) + dos(i+1)) * de / 2.d0
     if(i+1 > numdos) then
        write(*,'(a)') "ERROR: Energy scale in DOS does not cross zero!!"
        stop
     endif
  enddo
  write(*,'(a37,2x,f13.6)') "Integration of DOS up to Fermi level:", sum1
  dos = dos * nelect / sum1
  sum1 = 0.d0
  i=0
  do
     i = i+1
     if(edos(i+1) > 0.d0) exit
     de = (edos(i+1) - edos(i))
     sum1 = sum1 + (dos(i) + dos(i+1)) * de / 2.d0
     if(i+1 > numdos) then
        write(*,'(a)') "ERROR: Energy scale in DOS does not cross zero!!"
        stop
     endif
  enddo
  write(*,'(a37,2x,f13.6)') "Renormalised integrated DOS         :", sum1
  write(*,*)
!
! Now do algorithm to determine self-consistent Fermi level. Begin scanning
! at E_F = 0, go forwards by 1/10 of difference between E_F = 0 and E_F = 
! emax. Determine difference of (*) from zero, when sign of difference changes,
! reverse direction of scan reducing step size. Continue until step
! size is lower than a particular convergence criterion, or the difference
! is less than same. Take care for large exponentials (which return infinity),
! if they occur set the calculated value at a max. If this max is maintained
! between more than two steps, reverse the direction of the scan.
!
!
! Calculate neutral defect concentrations
!
  if(ndefects > 0) then
     conc0 = 0.d0
     do i=1,ndefects
        check0 = .true.
        do j=1,ncharge(i)
           if(abs(charge(i,j)) < small) then
              ener0(i) = energy(i,j)
              deg0(i) = deg(i,j)
              check0 = .false.
           endif
        enddo
        if(check0) cycle
        expfac = -ener0(i) / (kboltz * temperature)
        if(expfac > largefac) then
           conc0(i) = vlarge
        else
           conc0(i) = nsite(i) * deg0(i) * exp(expfac)
        endif
     enddo
  endif
!
! Set initial values for energy scan
!
  estep = emax / 10.d0
  direction = 1.d0
  delphi = -estep
!
  i = 0
  kmax = 0
  kmin = 0
  check_under = .false.
  check_under2 = .false.
  checkconv = .false.
!
  do
!
! Step E_F (don't step E_F if convergence in estep has been found in the
! previous iteration)
!
     if(checkconv .eqv. .false.) delphi = delphi + direction * estep
!
! Get n0 and p0 using integrals (equations 28.9 in Ashcroft Mermin)
!
     n0 = 0.d0
     p0 = 0.d0
     do j=1,numdos-1
        if(edos(j+1) <= 0.d0) then
           de = (edos(j+1) - edos(j))
           p0 = p0 + ( dos(j) / ( 1 + exp((delphi - edos(j)) / (kboltz * temperature)) )&
                & + dos(j+1) / (1 + exp((delphi - edos(j+1)) / (kboltz * temperature)) ) )&
                * de / 2.d0
        endif
        if(edos(j) >= egap) then
           de = (edos(j+1) - edos(j))
           n0 = n0 + ( dos(j) / ( 1 + exp((edos(j) - delphi) / (kboltz * temperature)) )&
                & + dos(j+1) / (1 + exp((edos(j+1) - delphi) / (kboltz * temperature)) ) )&
                * de / 2.d0
        endif
     enddo
     lhs = p0
     rhs = n0
!
! Get defect concentrations at E_F
!
     if(ndefects > 0) then
        do j=1,ndefects
!
! Get concentrations for each charge
!
           do k=1,ncharge(j)
              if(abs(charge(j,k)) < small) cycle 
              enerq = energy(j,k) + charge(j,k) * delphi
              expfac = -enerq / (kboltz * temperature)
              if(expfac > largefac) then
                 concd(j,k) = vlarge
              else
                 concd(j,k) = nsite(j) * deg(j,k) * exp(expfac)
              endif
              if(charge(j,k) < small) then
                 rhs = rhs + concd(j,k) * abs(charge(j,k))
              else
                 lhs = lhs + concd(j,k) * abs(charge(j,k))
              endif
           enddo
        enddo
     endif
!
! Get diff, check how small it is, get direction to go to reduce it
!
     diff = rhs - lhs
!
! If estep was converged in the previous iteration, exit now after recalculating
! all concentrations and diff (see 19/4/18)
!
     if(checkconv) exit
!
! Need to check for real underflow in exponentials, which is a problem if
! both the rhs and lhs are identically zero. If all terms become
! zero due to underflow, we find the range of fermi level within which
! this happens, and as an approximation return the midpoint as the solution
!
     if(rhs == 0.d0 .and. lhs == 0.d0) then
        if(check_under) cycle
        check_under = .true.
        write(*,'(a)') "WARNING: real number underflow due to v. small exponentials"
        write(*,'(a)') "         - treat results with caution!!"
        delphi_min = delphi
        if(estep * 1.d-4 > small) then
           estep = estep * 1.d-4
        else
           estep = small
        endif
     else
        if(check_under) then
           if(check_under2) then
!
! BUGFIX 23/4/18****
! After getting midpoint of range where everything went to zero as solution,
! need to go back and recalculate concentrations and diff. Otherwise, in very
! rare cases can get concentrations that are too large (= vlarge, but when 
! multiplied by 10^24 to get correct units, end up with 'infinity'), as the
! concentrations were calculated after exponentials were no longer zero
!
              delphi = 0.5d0 * (delphi - direction * estep + delphi_min)
              checkconv = .true.
              cycle
           else
              check_under2 = .true.
              temp1 = delphi
              delphi = delphi_min - direction * estep
              delphi_min = temp1
              direction = direction * (-1.d0)
           endif
        else
!
! Otherwise all is ok and we can perform the normal procedure
!
! BUGFIX APR 2018: put an else clause before this block, to 
! avoid the situation where you get a direction change directly 
! after a direction change when check_under2 is set to true
! (18/4/18)
!
           if(abs(diff) <= small*(rhs + lhs)) exit
!
! BUGFIX SEP 2018 (see 14-9-18): removed block of code that was originally
! put in place to check if diff = diffold and diff >= vlarge. The method was
! cumbersome, relying on many checks, and failed sometimes. Now the code 
! deals with vlarge diffs in the following way:
! 1. If |diff| > vlarge but |diffold| < vlarge, delphi goes back to previous
!    point and steps in same direction as before with smaller estep (to 
!    find where diff gets too big along a finer grid and change direction
!    appropriately)
! 2. If |diff| > |diffold| but both are > vlarge, direction is changed and
!    step size is reduced, even if it occurs in the first two steps (to avoid
!    going off into a region where diff > vlarge and getting stuck)
! 3. If |diff| < |diffold|, both are > vlarge, but there is a sign change
!    between them (checked by dividing one by the other and seeing if its < 0)
!    the direction is changed and step size reduced as the solution lies 
!    somewhere between them (the same situation but with |diff| > |diffold| is
!    covered in 2 above)
!
           if(i > 0) then 
              if(abs(diff) > abs(diffold)) then
                 if(abs(diff) >= vlarge .and. abs(diffold) < vlarge) then
                    delphi = delphi - direction * estep
                    diff = diffold
                    estep = estep / 4.d0
                 else
                    direction = direction * (-1.d0)
                    if(i > 1 .or. abs(diff) >= vlarge) estep = estep / 10.d0
                 endif
              else
                 if( (abs(diff) > vlarge .and. abs(diffold) > vlarge) .and. &
                      (diff / diffold < 0.d0) ) then
                    direction = direction * (-1.d0)
                    estep = estep / 10.d0
                 endif
              endif
           endif
        endif
        if(estep < small) then
!
! The solution is deemed converged. However, we must take a step back and
! recalculate the concentrations at the previous delphi, as these will be
! closer to the true solution (diff is lower in the previous step). Remember
! that direction has been changed, so we add rather than subtract direction *
! estep * 10.d0 (the 10 accounts for the decrease in step size between steps)
! See 19/4/18
!
           checkconv = .true.
           delphi = delphi + direction * estep * 10.d0 
           exit
        endif
        i = i+1
        diffold = diff
!
! If scan gets to boundaries of energy scale reverse direction, but only once
!
        if(delphi > emax) then
           if(kmax == 0) then
              direction = direction * (-1.d0)
              kmax = kmax+1
           else
              write(*,'(a)') "Sorry, no solution found within range of DOS :("
              stop
           endif
        endif
        if(delphi < emin) then
           if(kmin == 0) then
              direction = direction * (-1.d0)
              kmin = kmin+1
           else
              write(*,'(a)') "Sorry, no solution found within range of DOS :("
              stop
           endif
        endif
     endif
  enddo
!
! Write results - E_F and concentrations
!
  write(*,'(a)') "Solution found!"
  write(*,*)
  write(*,'(a41,x,e22.13e3,2x,a5)') "Condition (n + acceptors) - (p + donors):", &
       &diff * 1.d24 / volume, "cm^-3"
  write(*,*) 
  write(*,'(a)') "Results:"
  write(*,'(a)') "--------"
  write(*,'(a16,x,f21.13,2x,a4)') "SC Fermi level :", delphi, "(eV)"
  write(*,*)  
  write(*,'(a)') "Concentrations:"
  write(*,'(a16,x,e22.13e3,2x,a5)') "n (electrons)  :", n0 * 1.d24 / volume, "cm^-3" 
  write(*,'(a16,x,e22.13e3,2x,a5)') "p (holes)      :", p0 * 1.d24 / volume, "cm^-3" 
  if(ndefects > 0) then
     do i=1,ndefects
        concall(i) = conc0(i)
        do j=1,ncharge(i)
           concall(i) = concall(i) + concd(i,j)
        enddo
        write(*,'(a10,a6,x,e22.13e3,2x,a5)') name(i), "     :", concall(i) * 1.d24 / &
             &volume, "cm^-3"
     enddo
     write(*,*)
     write(*,'(a)') "Breakdown of concentrations for each defect charge state:"
     do i=1,ndefects
        write(*,'(a)') "--------------------------------------------------------"
        if(concall(i) == 0.d0) then
           write(*,'(a10,a6,3x,a34)') name(i), "     :", "Zero total - &
                &cannot give breakdown"
           cycle
        endif
        write(*,'(a10,a6,3x,a6,3x,a20,x,a7)') name(i), "     :", "Charge", &
             &"Concentration(cm^-3)", "% total"
        do j=1,ncharge(i)
           if(abs(charge(i,j)) < small) then
              write(*,'(15x,a1,3x,f5.1,3x,e22.13e3,x,f7.2)') ":", charge(i,j), conc0(i) * &
                   1.d24 / volume, conc0(i) * 100.d0 / concall(i)
           else
              write(*,'(15x,a1,3x,f5.1,3x,e22.13e3,x,f7.2)') ":", charge(i,j), concd(i,j) * &
                   1.d24 / volume, concd(i,j) * 100.d0 / concall(i)
           endif
        enddo
     enddo
  endif
!
end program sc_fermi
