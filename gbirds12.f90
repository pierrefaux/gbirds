program gbirds
!
!This program is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!any later version.

!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.

!You should have received a copy of the GNU General Public License
!along with this program.  If not, see <http://www.gnu.org/licenses/>.

!Pierre FAUX (pierrefaux@gmail.com), 2019
!Version 1.2, last update on 2019-10-15
!
implicit none
type i1cell
        integer:: i, n, p
        integer*1, allocatable:: v(:,:)
end type
type icell
        integer:: i, n, p
        integer, allocatable:: v(:)
end type
type icell2
        integer:: i, n, p
        integer, allocatable:: v(:,:)
end type
type sample
        integer:: i, u, n, nr, pop, nsg, nur, &
                bclength, bcvalue, lane, lowq
        logical:: kept
        integer, allocatable:: libp(:), libd(:), pos(:)
        integer*1, allocatable:: isu(:)
        character*10:: barcode
        character*20:: id
end type
integer:: h, i, j, k, kk, l, m, n, u, v, lg, lg1, nlib, nt, &
        nbgr, g, gr, nu, c, maxsize, ntag, khash, lhash, nbdif, nsplit, ntag0, &
        sp, ntagmax, ndi, nindel, nindel0, nbpb, npi, js, je0, je1, r, nsg, &
        trfn, nind, nlib1, ind, ind_nmax, lib_nmax, sg_nmax, t, nureads, &
        nthreads, usflag, tdepthmin, nureads0, sel, nel, npats, &
        nfemvec(256), nfem, nmal, best1, best2, maxnureads, GCvec(256), &
        best11, optim_opt, ntopocc, nureadsg, threstopocc, nupats, &
        lane, nlane, raw_nmax, fourtoten, lgtot, jj, lanelgflag, nlib0, &
        phredthr, phredthravg1, phredthravg2, modcheck, lp2, lp1, minureads, &
        naddureadsg, nind0, depthmin, covmed, covsin, khash2, lhash2, &
        lmark, np, nm1, nm2, nm3, nm4
integer, allocatable:: w2(:), w3(:,:), w4(:), p(:), w5(:,:), depth(:,:), &
        w6(:), w8(:), singletons(:), xref(:,:), tdepth(:), &
        nbocc(:), nindcov(:), topocc(:,:), sgref(:), indsM(:), indsF(:), &
        vcheck(:), icheck1(:), icheck2(:), icheck3(:), icheck4(:), qual(:), &
        samples(:), addsgref(:), whash2(:), libseqir(:,:), liborigin(:,:), &
        tags(:,:), indels(:,:), redunds(:,:), marker(:,:), pairs(:,:), &
        AC(:,:,:), AD(:,:,:), GC(:,:), GD(:,:), GCII(:,:), GCIII(:,:,:), &
        GCIV(:,:)
integer*1:: l1, ii
integer*1, allocatable:: seq(:), libseqi(:,:), w1i1(:,:), &
        compmat(:,:,:), ureads(:,:), ureads1(:,:), w2i1(:), &
        w3i1(:,:), patterns(:,:), upatterns(:,:), sexv(:), ureadsg(:,:), &
        raws(:,:), GCout(:), allraws(:,:), addureadsg(:,:), seqivec(:,:)
integer*2, allocatable:: bcid(:,:), indraws(:)
integer*8:: nbc, kkhash2
logical:: save_opt, interspec_check
logical, allocatable:: isW(:)
real:: t1, t0, avg, std, nsexminpc, maxsc1, maxsc2, maxsc11, maxsc21, bounds(2)
real, allocatable:: rw2(:), scores(:,:), cpm(:,:)
real*8:: sexrate, bestsexrate
real*8, allocatable:: lpt(:)
character*1:: re(4), barcode(10), junk1
character*1, allocatable:: seq0(:), seq1(:), phred0(:), seqirvec(:,:)
character*3:: current_version
character*6:: sexc(2)
character*10:: junk10
character*20:: fmt1, fmt2, junk20
character*50:: popfile, listfile, parfile, junk50
character*100:: lanefile
character*200:: junk200
character*1000:: junk1000
type(icell), allocatable:: newg(:), oldg(:,:), split(:,:), lanes(:)
type(icell2), allocatable:: recip(:)
type(sample), allocatable:: indiv(:)

current_version='1.2'

call setup_parameters
call setup_GCvec
call setup_nfemvec
call setup_translation_vec

! System parameters
re(1)='T'; re(2)='G'; re(3)='C'; re(4)='A'
fourtoten=1048576+1
khash=4
npi=3
raw_nmax=20E6
ind_nmax=1E6
lib_nmax=6E6
sg_nmax=20E6
nsexminpc=15.0
ntagmax=5E5
khash2=11
kkhash2=4**(khash2+2)
allocate(whash2(khash2))
do i=1, khash2
        whash2(i)=4**(i-1)
enddo


! STEP 1: LANE PARSING

! STEP 1.1 : READING THE POPULATION MAP
junk50='LOADING OF POPULATION MAP'
call print_section_title(junk50)
call cpu_time(t0)
call parse_population_map
call cpu_time(t1)
print *, ' population map loaded, ', t1-t0, ' seconds '
if (usflag==0) then
        ! STEP 1.2 : SEQUENTIAL READING OF THE SEQUENCING LANES
        junk50='SEQUENCING LANES PARSING'
        call print_section_title(junk50)
        call cpu_time(t0)
        call parse_sequencing_lanes
        call cpu_time(t1)
        print *, '  all sequencing lanes parsed, elapsed time = ', t1-t0
        ! STEP 2 : UNIQUE SORT OF READS
        
        ! STEP 2.1 : PER SAMPLE
        junk50='UNIQUE SORT OF READS'
        call print_section_title(junk50)
        call cpu_time(t0)
        call sort_raws_per_sample
        call cpu_time(t1)
        print *, '  reads uniquely sorted for all samples, elapsed time = ', t1-t0
        
        ! STEP 2.2 : UNIQUE SORT OF POOLED READS 
        
        call cpu_time(t0)
        call sort_all_unique_raws
        call cpu_time(t1)
        print *, '  unique reads from samples pooled together and uniquely sorted, '
        print *, '  elapsed time = ', t1-t0 
        
        ! STEP 3 : QUALITY CONTROLS (READS AND SAMPLES)
        junk50='QUALITY CONTROL'
        call print_section_title(junk50)
        call cpu_time(t0)
        call quality_check
        call cpu_time(t1)
        print *, '  quality checks on reads and samples, '
        print *, '  elapsed time = ', t1-t0

        ! STEP 4 : WRITING .us AND .usg files
        junk50='WRITING BINARY FILES'
        call print_section_title(junk50)
        call cpu_time(t0)
        call save_us_files
        call cpu_time(t1)
        print *, '  unique reads and associated depths saved to .us and .usg files, '
        print *, '  elapsed time = ', t1-t0
else
        ! STEP 2b : LOADING OF .us AND .usg files
        junk50='LOADING BINARY FILES'
        call print_section_title(junk50)
        call cpu_time(t0)
        call load_binary_files
        call cpu_time(t1)
        print *, '  all binary files load in memory, '
        print *, '  elapsed time = ', t1-t0
endif

! STEP 4 : SEX IDENTIFICATION
junk50='SEX IDENTIFICATION'
call print_section_title(junk50)
call cpu_time(t0)
call apply_depth_filters
call sex_identification
if (save_opt) call save_library
call cpu_time(t1)
print *, '  sex identification of samples achieved, '
print *, '  elapsed time = ', t1-t0

! STEP 5 : CLUSTERING OF READS
junk50='CLUSTERING OF SEQUENCES'
call print_section_title(junk50)
call cpu_time(t0)
call get_libseqi        ! ureads must be un-hashed for marker search
call compar_reads
call cpu_time(t1)
print *, '  clustering of reads into homolog sequences achieved, '
print *, '  elapsed time = ', t1-t0

! STEP 6 : GENOTYPE CALLING
junk50='GENOTYPE CALLING'
call print_section_title(junk50)
call cpu_time(t0)
call biAssemble2
call cpu_time(t1)
print *, '  calling SNPs and indels from sequence achieved, '
print *, '  elapsed time = ', t1-t0

! STEP 7 : GENOTYPE FILTERING
junk50='GENOTYPE FILTERING'
call print_section_title(junk50)
call cpu_time(t0)
call filtering
call cpu_time(t1)
print *, '  filtering genotypes, '
print *, '  elapsed time = ', t1-t0

contains

subroutine parse_population_map
implicit none
open(newunit=u, file=trim(adjustl(listfile)))
nlane=count_lines(u)
print *, 'number of sequencing lanes = ', nlane
allocate(lanes(nlane))
close(u)
! Initialization of the indiv structure
open(newunit=u, file=trim(adjustl(popfile)))
nind=count_lines(u)
print *, ' #indiv = ', nind
nind0=nind
allocate(indiv(nind),bcid(fourtoten,nlane))
bcid=0
do i=1, nind
        indiv(i)%i=i
        read(u,*) junk10, junk20, indiv(i)%pop, indiv(i)%lane
        j=index(junk10,' ')-1
        indiv(i)%bclength=j
        read(junk10(:j),'(a)') indiv(i)%barcode
        barcode(:)='A'
        do k=1, j
                barcode(k)=junk10(k:k)
        enddo
        do k=j+1, j+4
                if (k>10) exit
                barcode(k)=re(k-j)
        enddo
        call setup_bcid(i,j,barcode)
        j=index(junk20,' ')-1
        indiv(i)%id=junk20(:j)
        indiv(i)%kept=.true.
        indiv(i)%lowq=0
        print *, i, indiv(i)%barcode, indiv(i)%bclength, indiv(i)%bcvalue, &
                indiv(i)%lane
enddo
close(u)
end subroutine parse_population_map

subroutine parse_sequencing_lanes
implicit none
real:: t0, t1
call cpu_time(t1)
open(newunit=u, file=trim(adjustl(listfile)))
do lane=1, nlane
        if (lanelgflag==0) then
                read(u,*) lanefile
        else
                read(u,*) lanefile, i
                lanes(lane)%p=i
        endif
        open(newunit=lanes(lane)%i, file=trim(adjustl(lanefile)))
        lanes(lane)%n=count(indiv(:)%lane==lane)
        allocate(lanes(lane)%v(lanes(lane)%n))
        k=0
        do i=1, nind
                if (indiv(i)%lane/=lane) cycle
                k=k+1
                lanes(lane)%v(k)=i
        enddo
enddo
rewind(u)
! Determination of the length of the reads
read(lanes(1)%i,*) 
read(lanes(1)%i,'(a1000)') junk1000
rewind(lanes(1)%i)
lgtot=index(junk1000,' ')-1
print *, '       raw read length = ', lgtot
lg=lgtot-(maxval(indiv(:)%bclength)+4) !! this is assuming that the RE lenght is always 4
print *, '      good read length = ', lg
lhash=ceiling(real(lg)/real(khash))
print *, '    hashed read length = ', lhash
lg1=lhash*khash
print *, ' un-hashed read length = ', lg1
lp1=lg+khash
lp2=lg1+khash
! re-open with recl=lgtot (for speed consideration) and count the number of reads
nlib0=0
do lane=1, nlane
        close(lanes(lane)%i)
        read(u,*) lanefile
        open(newunit=lanes(lane)%i, file=trim(adjustl(lanefile)), recl=lgtot)
        if (lanelgflag==0) then
                lanes(lane)%p=count_lines(lanes(lane)%i)/4
        endif
        nlib0=nlib0+lanes(lane)%p
enddo
print *, ' total number of reads in lanes = ', nlib0
call cpu_time(t0)
print *, ' elapsed time for obtaining size of sequencing lanes, ', &
         t0-t1, ' seconds'

! Reading/parsing the sequencing lanes
call setup_phred_check
allocate(raws(nlib0,lhash), seq0(lgtot), seq1(lg1), &
        indraws(nlib0), phred0(lgtot), qual(lp2))
!indraws=0 ! in order to reduce memory pressure, we do not initialize this
do i=1, lg1
        seq1(i)='0'
enddo
qual=0
write(fmt1,'(a,i5,a)') '(',lgtot,'a1)'
! the following loop should be parallelizable
nlib=0
do lane=1, nlane
        call cpu_time(t0)
        nlib1=nlib+1
        do i=1, lanes(lane)%p
                read(lanes(lane)%i,*)                 ! part1
                read(lanes(lane)%i,fmt1) seq0        ! part2
                ! 1st check: does the barcode part contain Ns?
                kk=count(seq0(1:10)=='N')
                if (kk>0) then
                        read(lanes(lane)%i,*)   ! part 3
                        read(lanes(lane)%i,*)   ! part 4
                        cycle
                endif
                ! 2nd check: can we identify the sample?
                j=translate_barcode_10(seq0(1:10))
                k=bcid(j,lane)
                if (k==0) then
                        read(lanes(lane)%i,*)   ! part 3
                        read(lanes(lane)%i,*)   ! part 4
                        cycle
                endif
                ! 3rd check: does the read (incl. RE) reach required quality?
                read(lanes(lane)%i,*)         ! part3
                read(lanes(lane)%i,fmt1) phred0        ! part4
                qual(1:lp1)=iachar( &
                        phred0(indiv(k)%bclength+1:indiv(k)%bclength+4+lg))
                if (phred_check(qual,lp2)) then
                        indiv(k)%lowq=indiv(k)%lowq+1
                        cycle
                endif
                ! if all check are OK, keep the read
                nlib=nlib+1
                indraws(nlib)=k
                seq1(1:lg)=seq0(indiv(k)%bclength+5:indiv(k)%bclength+4+lg)
                !! this is assuming that the RE lenght is always 4
                raws(nlib,:)=translate_reads_3(seq1,lg1,lhash)
        enddo
        call cpu_time(t1)
        print *, ' lane # ', lane, ' processed, ', t1-t0, ' seconds'
        print *, '   # of reads assigned in this lane = ', nlib-nlib1+1
        print *, '    unassigned reads # in this lane = ', &
                lanes(lane)%p-(nlib-nlib1+1)
        print *, '       SAMPLE_ID       #GOOD_READS       #BAD_READS'
        allocate(w2(nlib))
        w2=(/(i, i=1, nlib)/)
        do i=1, lanes(lane)%n
                j=lanes(lane)%v(i)
                indiv(j)%nr=count(indraws(1:nlib)==j)
                allocate(indiv(j)%pos(indiv(j)%nr))
                indiv(j)%pos=pack(w2,indraws(1:nlib)==j)
                print *, '      ', indiv(j)%id, indiv(j)%nr, &
                        indiv(j)%lowq
        enddo
        deallocate(w2)
enddo
print *, ' actual number of reads kept for next steps = ', nlib
deallocate(seq0)
end subroutine parse_sequencing_lanes

subroutine setup_parameters
implicit none
integer:: i, j, k, narg, modarg, narg2, nopt
logical:: flagopt
character*50:: argname, argvalue
! List of default values for options
        sexc(1)='male'
        sexc(2)='female'
        minureads=5E5
        usflag=1
        phredthr=40
        lanelgflag=1
        save_opt=.false.
        interspec_check=.true.
        optim_opt=0
        sexrate=.5
        covmed=6
        covsin=10
        tdepthmin=3
        depthmin=3
        nbdif=2
        nbpb=9
narg=command_argument_count()
modarg=mod(narg,2)
if (narg<1.or.modarg/=0) then
        print *, ' ERROR :: something wrong in argument list '
        stop
endif
if (narg<4) then
        print *, ' ERROR :: mandatory arguments missing '
        stop
endif
narg2=narg/2
flagopt=.false.
do i=1, narg2
        k=i*2
        j=k-1
        call getarg(j,argname)
        call getarg(k,argvalue)
        select case (trim(adjustl(argname)))
        case("--map")
                popfile=argvalue
        case("--files")
                listfile=argvalue
        case("--options")
                flagopt=.true.
                parfile=argvalue
        case default
        endselect
enddo
call print_section_separation
print *, '                                          ****** '
print *, '                                     GBIRDS version ',current_version
print *, '                                  (Faux and Santos, 2019)'
print *, '                                          ****** '
print *, '                             contact us: pierrefaux@gmail.com'
print *, '                                          ****** '
print *, ''
print '(a,a)', '        Population map in ', trim(adjustl(popfile))
print '(a,a)', '        List of FASTQ files in ', trim(adjustl(listfile))
if (flagopt) then
print '(a,a)', '        Options passed in file ', trim(adjustl(parfile))
open(newunit=u, file=trim(adjustl(parfile)))
nopt=count_lines(u)
do i=1, nopt
        read(u,*) argname, argvalue
        select case (trim(adjustl(argname)))
        case("qual_threshold")
                argvalue=trim(adjustl(argvalue))
                read(argvalue,'(i2)') phredthr
                print '(a,a,i3)', '                        ', &
                        '* BP quality threshold set to ', phredthr
        case("min_sample_coverage")
                argvalue=trim(adjustl(argvalue))
                read(argvalue,'(i4)') depthmin
                print '(a,a,i3)', '                        ', &
                        '* Minimum sequence coverage per sample set to ', depthmin
        case("min_total_coverage")
                argvalue=trim(adjustl(argvalue))
                read(argvalue,'(i4)') tdepthmin
                print '(a,a,i3)', '                        ', &
                        '* Minimum total sequence coverage set to ', tdepthmin
        case("min_median_coverage")
                argvalue=trim(adjustl(argvalue))
                read(argvalue,'(i4)') covmed
                print '(a,a,i3)', '                        ', &
                        '* Minimum median coverage of alleles set to ', covmed
        case("min_singleton_allele_coverage")
                argvalue=trim(adjustl(argvalue))
                read(argvalue,'(i2)') covsin
                print '(a,a,i3)', '                        ', &
                        '* Minimum coverage of singleton allele set to ', covsin
        case("start_from_us")
                argvalue=trim(adjustl(argvalue))
                read(argvalue,'(i1)') usflag
                if (usflag<0.or.usflag>1) then
                        print '(a,a,i1)', '                        ', &
                                ' ERROR: value out of range : ', usflag
                        stop
                endif
                if (usflag==1) then
                        print '(a,a)', '                        ', & 
                                '* Analysis will start from existing *.us'
                else
                        print '(a,a)', '                        ', &
                                '* Analysis will start from fastq files'
                endif
        case("lane_lengths_given")
                argvalue=trim(adjustl(argvalue))
                read(argvalue,'(i1)') lanelgflag
                if (lanelgflag<0.or.lanelgflag>1) then
                        print '(a,a,i1)', '                        ', &
                                ' ERROR: value out of range : ', lanelgflag
                        stop
                endif
                if (lanelgflag==1) then
                        print '(a,a)', '                        ', & 
                                '* Number of reads in fastq passed in file list'
                else
                        print '(a,a)', '                        ', &
                                '* Reads in fastq have to be numbered'
                endif
        case("min_ureads")
                argvalue=trim(adjustl(argvalue))
                read(argvalue,*) minureads
                if (minureads<=0) then
                        print '(a,a,i9)', '                        ', &
                                ' WARNING : value out of range : ', minureads
                endif
                print '(a,a,i9)', '                        ', &
                        '* Minimum nb. of reads set to ', minureads
        case("save_library")
                argvalue=trim(adjustl(argvalue))
                read(argvalue,*) j
                if (j==1) then
                        save_opt=.true.
                        print '(a,a)', '                        ', &
                                '* Save reads library to files'
                elseif (j==0) then
                        save_opt=.false.
                else
                        save_opt=.false.
                        print '(a,a,i9)', '                        ', &
                                ' WARNING : value out of range : ', j
                endif
        case("interspecific_check")
                argvalue=trim(adjustl(argvalue))
                read(argvalue,*) j
                if (j==0) then
                        interspec_check=.false.
                        print '(a,a)', '                        ', &
                                '* No interspecific samples check'
                elseif (j==1) then
                        interspec_check=.true.
                else
                        interspec_check=.false.
                        print '(a,a,i9)', '                        ', &
                                ' WARNING : value out of range : ', j
                endif
        case("optimize_prior_rate")
                argvalue=trim(adjustl(argvalue))
                read(argvalue,*) optim_opt
                if (optim_opt<0.or.optim_opt>1) then
                        print '(a,a,i9)', '                        ', &
                                ' WARNING : value out of range : ', optim_opt
                endif
                print '(a,a)', '                        ', &
                        '* Optimize prior sex distribution rate'
        case("prior_sex_rate")
                argvalue=trim(adjustl(argvalue))
                read(argvalue,*) j
                sexrate=real(j)/100.d0
                if (j<0.or.j>100) then
                        print '(a,a,i1)', '                        ', &
                                ' ERROR : value out of range for prior_sex_rate : ', j
                endif
                print '(a,a,i3)', '                        ', &
                        '* Prior sex distribution rate set to ', j
        case("max_nb_snps")
                argvalue=trim(adjustl(argvalue))
                read(argvalue,*) j
                if (j<0.or.j>3) then
                        print '(a,a,i5)', '                        ', &
                                ' ERROR : value out of range for max_nb_snps : ', j
                        stop
                endif
                print '(a,a,i1)', '                        ', &
                        '* Maximum number of SNPs per read set to ', j
                nbdif=j
        case("max_nb_indels")
                argvalue=trim(adjustl(argvalue))
                read(argvalue,*) j
                if (j<0.or.j>9) then
                        print '(a,a,i5)', '                        ', &
                                ' ERROR : value out of range for max_nb_indels : ', j
                        stop
                endif
                print '(a,a,i1)', '                        ', &
                        '* Maximum number of indels per read set to ', j
                nbpb=j
        case default
                print '(a,a,a)', '                        ', &
                        '* (Unrecognized option: ', trim(adjustl(argvalue))
        endselect
enddo
else
print '(a,a)', '        (No extra options passed) '
endif
print *, ''
end subroutine setup_parameters

function translate_reads_3(s,l,lh) result(v)
! the length lhash is a multiple of khash !!
! this version is ACTG: A=0, C=1, G=3, T=2
implicit none
integer:: i, j, k, l, lh
integer*1:: v(lh)
character*1:: s(l)
v=0
k=0
do i=1, lhash
        do j=1, 3
                k=k+1
                v(i)=v(i)+ibits(iachar(s(k)),1,2)
                v(i)=ishft(v(i),2)
        enddo
        k=k+1
        v(i)=v(i)+ibits(iachar(s(k)),1,2)
enddo
end function translate_reads_3

function translate_barcode_10(bc) result(v)
! this version is ACTG: A=0, C=1, G=3, T=2
implicit none
integer:: j, v
character*1:: bc(10)
v=0
do j=1, 9
        v=v+ibits(iachar(bc(j)),1,2)
        v=ishft(v,2)
enddo
v=v+ibits(iachar(bc(10)),1,2)+1
end function translate_barcode_10


function count_lines(un_file) result(n)
implicit none
integer:: un_file, n, k
character*1:: junk
rewind(un_file); n=0
do
        read(un_file,'(a1)',iostat=k)junk
        if(k/=0)exit
        n=n+1
enddo
rewind(un_file)
end function count_lines

subroutine setup_bcid(ind,l,barcode)
implicit none
integer:: i, j, k, ind, l
character*1:: barcode(10), alphabet(4), barcode2(10)
alphabet(1)='A'
alphabet(2)='C'
alphabet(3)='T'
alphabet(4)='G'
! the following has the form barcode + RE + trailing A's
! it is the canonical value for barcode
indiv(ind)%bcvalue=translate_barcode_10(barcode)
! the next ones are the possible variations with one bp error (4*length of
! barcode)
do i=1, l
        do k=1, 4
                barcode2=barcode
                barcode2(i)=alphabet(k)
                j=translate_barcode_10(barcode2)
                bcid(j,indiv(ind)%lane)=ind
        enddo
enddo
! the next ones are only for 5bp-long barcode
if (l==5) then
        do k=2, 4
                barcode2=barcode
                barcode2(10)=alphabet(k)
                j=translate_barcode_10(barcode2)
                bcid(j,indiv(ind)%lane)=ind
        enddo
endif
! for security, the canonical value is set after others in bcid
bcid(indiv(ind)%bcvalue,indiv(ind)%lane)=ind
end subroutine setup_bcid

subroutine setup_phred_check
implicit none
integer:: i, j, k, d
d=lhash+1
modcheck=mod(lg,khash)
allocate(vcheck(d), icheck1(d), icheck2(d), &
        icheck3(d), icheck4(d))
phredthravg1=khash*phredthr
if (modcheck==0) then
        phredthravg2=phredthravg1
else
        phredthravg2=modcheck*phredthr
endif
j=(d-1)*4
k=1; icheck1=(/(i, i=k,j+k,4)/)
k=2; icheck2=(/(i, i=k,j+k,4)/)
k=3; icheck3=(/(i, i=k,j+k,4)/)
k=4; icheck4=(/(i, i=k,j+k,4)/)
end subroutine setup_phred_check

function phred_check(phredseq,l) result (f)
! f is a flag
!       = 1 if the read does not reach the required quality
!       = 0 otherwise (good read)
implicit none
integer:: l, phredseq(l)
logical:: f
vcheck=phredseq(icheck1)+phredseq(icheck2)+ &
        phredseq(icheck3)+phredseq(icheck4)
f=minval(vcheck(1:lhash))<phredthravg1.or.vcheck(lhash+1)<phredthravg2
end function phred_check

subroutine sort_raws_per_sample
implicit none
real:: t0, t1
integer*1, allocatable:: libseqir1(:,:), db(:)

nt=1
allocate(oldg(ind_nmax,nt), singletons(sg_nmax))
do ind=1, nind
        call cpu_time(t0)
        t=1
        nlib1=indiv(ind)%nr
        allocate(db(nlib1), libseqir1(nlib1,lhash))
        libseqir1=raws(indiv(ind)%pos,:)
        js=1
        je0=1
        r=0
        nsg=0
        allocate(oldg(1,t)%v(nlib1))
        oldg(:,t)%n=0
        oldg(1,t)%v=(/(i,i=1,nlib1)/)
        oldg(1,t)%n=nlib1
        do i=1,lhash
                if (i>=npi) then ! work with kmer=2 !! optimum should be 2 or 3, i.e. it
                        do ii=4, 0, -4
                                db=ibits(libseqir1(:,i),ii,4)
                                je1=je0
                                if (js<=je0) then
                                        do j=js,je0
                                                if (oldg(j,t)%n>1) then
                                                        do l1=0,15
                                                                m=count(db(oldg(j,t)%v)==l1)
                                                                if (m>0) then
                                                                        je1=je1+1
                                                                        if (je1>ind_nmax) je1=1
                                                                        allocate(oldg(je1,t)%v(m))
                                                                        oldg(je1,t)%v=pack(oldg(j,t)%v,db(oldg(j,t)%v)==l1)
                                                                        oldg(je1,t)%n=m
                                                                endif
                                                        enddo
                                                else
                                                        nsg=nsg+1
                                                        singletons(nsg)=oldg(j,t)%v(1)
                                                endif
                                                deallocate(oldg(j,t)%v)
                                                oldg(j,t)%n=0
                                        enddo
                                else
                                r=r+1
                                        do j=js,ind_nmax
                                                if (oldg(j,t)%n>1) then
                                                        do l1=0,15
                                                                m=count(db(oldg(j,t)%v)==l1)
                                                                if (m>0) then
                                                                        je1=je1+1
                                                                        if (je1>ind_nmax) je1=1
                                                                        allocate(oldg(je1,t)%v(m))
                                                                        oldg(je1,t)%v=pack(oldg(j,t)%v,db(oldg(j,t)%v)==l1)
                                                                        oldg(je1,t)%n=m
                                                                endif
                                                        enddo
                                                else
                                                        nsg=nsg+1
                                                        singletons(nsg)=oldg(j,t)%v(1)
                                                endif
                                                deallocate(oldg(j,t)%v)
                                                oldg(j,t)%n=0
                                        enddo
                                        do j=1,je0
                                                if (oldg(j,t)%n>1) then
                                                        do l1=0,15
                                                                m=count(db(oldg(j,t)%v)==l1)
                                                                if (m>0) then
                                                                        je1=je1+1
                                                                        if (je1>ind_nmax) je1=1
                                                                        allocate(oldg(je1,t)%v(m))
                                                                        oldg(je1,t)%v=pack(oldg(j,t)%v,db(oldg(j,t)%v)==l1)
                                                                        oldg(je1,t)%n=m
                                                                endif
                                                        enddo
                                                else
                                                        nsg=nsg+1
                                                        singletons(nsg)=oldg(j,t)%v(1)
                                                endif
                                                deallocate(oldg(j,t)%v)
                                                oldg(j,t)%n=0
                                        enddo
                                endif
                                if (je0<ind_nmax) then
                                        js=je0+1
                                else        ! this case is very unlikely
                                        js=1
                                endif
                                je0=je1
                        enddo
                else ! work with kmer=1
                        do ii=6, 0, -2
                                db=ibits(libseqir1(:,i),ii,2)
                                je1=je0
                                if (js<=je0) then
                                        do j=js,je0
                                                if (oldg(j,t)%n>1) then
                                                        do l1=0,3
                                                                m=count(db(oldg(j,t)%v)==l1)
                                                                if (m>0) then
                                                                        je1=je1+1
                                                                        if (je1>ind_nmax) je1=1
                                                                        allocate(oldg(je1,t)%v(m))
                                                                        oldg(je1,t)%v=pack(oldg(j,t)%v,db(oldg(j,t)%v)==l1)
                                                                        oldg(je1,t)%n=m
                                                                endif
                                                        enddo
                                                else
                                                        nsg=nsg+1
                                                        singletons(nsg)=oldg(j,t)%v(1)
                                                endif
                                                deallocate(oldg(j,t)%v)
                                                oldg(j,t)%n=0
                                        enddo
                                else
                                r=r+1
                                        do j=js,ind_nmax
                                                if (oldg(j,t)%n>1) then
                                                        do l1=0,3
                                                                m=count(db(oldg(j,t)%v)==l1)
                                                                if (m>0) then
                                                                        je1=je1+1
                                                                        if (je1>ind_nmax) je1=1
                                                                        allocate(oldg(je1,t)%v(m))
                                                                        oldg(je1,t)%v=pack(oldg(j,t)%v,db(oldg(j,t)%v)==l1)
                                                                        oldg(je1,t)%n=m
                                                                endif
                                                        enddo
                                                else
                                                        nsg=nsg+1
                                                        singletons(nsg)=oldg(j,t)%v(1)
                                                endif
                                                deallocate(oldg(j,t)%v)
                                                oldg(j,t)%n=0
                                        enddo
                                        do j=1,je0
                                                if (oldg(j,t)%n>1) then
                                                        do l1=0,3
                                                                m=count(db(oldg(j,t)%v)==l1)
                                                                if (m>0) then
                                                                        je1=je1+1
                                                                        if (je1>ind_nmax) je1=1
                                                                        allocate(oldg(je1,t)%v(m))
                                                                        oldg(je1,t)%v=pack(oldg(j,t)%v,db(oldg(j,t)%v)==l1)
                                                                        oldg(je1,t)%n=m
                                                                endif
                                                        enddo
                                                else
                                                        nsg=nsg+1
                                                        singletons(nsg)=oldg(j,t)%v(1)
                                                endif
                                                deallocate(oldg(j,t)%v)
                                                oldg(j,t)%n=0
                                        enddo
                                endif
                                if (je0<ind_nmax) then
                                        js=je0+1
                                else        ! this case is very unlikely
                                        js=1
                                endif
                                je0=je1
                        enddo
                endif
        enddo
        if (js<=je0) then ! last scan to collect the last singletons
                do j=js,je0
                        if (oldg(j,t)%n==1) then
                                nsg=nsg+1
                                singletons(nsg)=oldg(j,t)%v(1)
                                deallocate(oldg(j,t)%v)
                                oldg(j,t)%n=0
                        endif
                enddo
        else
                do j=js,ind_nmax
                        if (oldg(j,t)%n==1) then
                                nsg=nsg+1
                                singletons(nsg)=oldg(j,t)%v(1)
                                deallocate(oldg(j,t)%v)
                                oldg(j,t)%n=0
                        endif
                enddo
                do j=1, je0
                        if (oldg(j,t)%n==1) then
                                nsg=nsg+1
                                singletons(nsg)=oldg(j,t)%v(1)
                                deallocate(oldg(j,t)%v)
                                oldg(j,t)%n=0
                        endif
                enddo
        endif
        indiv(ind)%n=count(oldg(:,t)%n>1)+nsg
!        print *, 'nbgr=', indiv(ind)%n
        allocate(indiv(ind)%libp(indiv(ind)%n),indiv(ind)%libd(indiv(ind)%n))
        indiv(ind)%nsg=nsg
        indiv(ind)%libp(1:nsg)=indiv(ind)%pos(singletons(1:nsg))
        indiv(ind)%libd(1:nsg)=1
        k=nsg
        do j=1, ind_nmax
                if (oldg(j,t)%n>=2) then
                        k=k+1
                        indiv(ind)%libp(k)=indiv(ind)%pos(oldg(j,t)%v(1))
                        indiv(ind)%libd(k)=oldg(j,t)%n
                endif
                if (allocated(oldg(j,t)%v)) deallocate(oldg(j,t)%v)
        enddo
        deallocate(libseqir1, db)
        maxsize=maxval(indiv(ind)%libd)
        call cpu_time(t1)
        print *, '    ', indiv(ind)%id(:len_trim(indiv(ind)%id)), &
                ' :: ', indiv(ind)%n, ' unique reads, sorted in ', &
                t1-t0, ' seconds'
enddo
deallocate(oldg, singletons)
end subroutine sort_raws_per_sample

subroutine sort_all_unique_raws
implicit none
integer:: i, j, k
integer*1, allocatable:: libseqir1(:,:), db(:)
real:: t0, t1
call cpu_time(t0)
nlib1=0
do i=1, nind
        if (indiv(i)%kept) nlib1=nlib1+indiv(i)%n
enddo
print *, ' # of individual-unique reads =', nlib1
allocate(libseqir1(nlib1,lhash), xref(nlib1,2), db(nlib1))
k=0
do i=1, nind
        if (.not.indiv(i)%kept) cycle
        j=indiv(i)%n
        libseqir1(k+1:k+j,:)=raws(indiv(i)%libp,:)
        xref(k+1:k+j,1)=i
        xref(k+1:k+j,2)=indiv(i)%libd
        indiv(i)%n=0
        deallocate(indiv(i)%libd, indiv(i)%libp)
        k=k+j
enddo
deallocate(raws,indraws)
call cpu_time(t1)
print *, ' elapsed time for setting xref=', t1-t0
allocate(oldg(lib_nmax,1), singletons(sg_nmax))
t=1
js=1
je0=1
r=0
nsg=0
allocate(oldg(1,t)%v(nlib1))
oldg(:,t)%n=0
oldg(1,t)%v=(/(i,i=1,nlib1)/)
oldg(1,t)%n=nlib1
do i=1,lhash
        if (i>=npi) then ! work with kmer=2 !! optimum should be 2 or 3, i.e. it
                do ii=4, 0, -4
                        db=ibits(libseqir1(:,i),ii,4)
                        je1=je0
                        if (js<=je0) then
                                do j=js,je0
                                        if (oldg(j,t)%n>1) then
                                                do l1=0,15
                                                        m=count(db(oldg(j,t)%v)==l1)
                                                        if (m>0) then
                                                                je1=je1+1
                                                                if (je1>lib_nmax) je1=1
                                                                allocate(oldg(je1,t)%v(m))
                                                                oldg(je1,t)%v=pack(oldg(j,t)%v,db(oldg(j,t)%v)==l1)
                                                                oldg(je1,t)%n=m
                                                        endif
                                                enddo
                                        else
                                                nsg=nsg+1
                                                singletons(nsg)=oldg(j,t)%v(1)
                                        endif
                                        deallocate(oldg(j,t)%v)
                                        oldg(j,t)%n=0
                                enddo
                        else
                        r=r+1
                                do j=js,lib_nmax
                                        if (oldg(j,t)%n>1) then
                                                do l1=0,15
                                                        m=count(db(oldg(j,t)%v)==l1)
                                                        if (m>0) then
                                                                je1=je1+1
                                                                if (je1>lib_nmax) je1=1
                                                                allocate(oldg(je1,t)%v(m))
                                                                oldg(je1,t)%v=pack(oldg(j,t)%v,db(oldg(j,t)%v)==l1)
                                                                oldg(je1,t)%n=m
                                                        endif
                                                enddo
                                        else
                                                nsg=nsg+1
                                                singletons(nsg)=oldg(j,t)%v(1)
                                        endif
                                        deallocate(oldg(j,t)%v)
                                        oldg(j,t)%n=0
                                enddo
                                do j=1,je0
                                        if (oldg(j,t)%n>1) then
                                                do l1=0,15
                                                        m=count(db(oldg(j,t)%v)==l1)
                                                        if (m>0) then
                                                                je1=je1+1
                                                                if (je1>lib_nmax) je1=1
                                                                allocate(oldg(je1,t)%v(m))
                                                                oldg(je1,t)%v=pack(oldg(j,t)%v,db(oldg(j,t)%v)==l1)
                                                                oldg(je1,t)%n=m
                                                        endif
                                                enddo
                                        else
                                                nsg=nsg+1
                                                singletons(nsg)=oldg(j,t)%v(1)
                                        endif
                                        deallocate(oldg(j,t)%v)
                                        oldg(j,t)%n=0
                                enddo
                        endif
                        if (je0<lib_nmax) then
                                js=je0+1
                        else        ! this case is very unlikely
                                js=1
                        endif
                        je0=je1
                enddo
        else ! work with kmer=1
                do ii=6, 0, -2
                        db=ibits(libseqir1(:,i),ii,2)
                        je1=je0
                        if (js<=je0) then
                                do j=js,je0
                                        if (oldg(j,t)%n>1) then
                                                do l1=0,3
                                                        m=count(db(oldg(j,t)%v)==l1)
                                                        if (m>0) then
                                                                je1=je1+1
                                                                if (je1>lib_nmax) je1=1
                                                                allocate(oldg(je1,t)%v(m))
                                                                oldg(je1,t)%v=pack(oldg(j,t)%v,db(oldg(j,t)%v)==l1)
                                                                oldg(je1,t)%n=m
                                                        endif
                                                enddo
                                        else
                                                nsg=nsg+1
                                                singletons(nsg)=oldg(j,t)%v(1)
                                        endif
                                        deallocate(oldg(j,t)%v)
                                        oldg(j,t)%n=0
                                enddo
                        else
                        r=r+1
                                do j=js,lib_nmax
                                        if (oldg(j,t)%n>1) then
                                                do l1=0,3
                                                        m=count(db(oldg(j,t)%v)==l1)
                                                        if (m>0) then
                                                                je1=je1+1
                                                                if (je1>lib_nmax) je1=1
                                                                allocate(oldg(je1,t)%v(m))
                                                                oldg(je1,t)%v=pack(oldg(j,t)%v,db(oldg(j,t)%v)==l1)
                                                                oldg(je1,t)%n=m
                                                        endif
                                                enddo
                                        else
                                                nsg=nsg+1
                                                singletons(nsg)=oldg(j,t)%v(1)
                                        endif
                                        deallocate(oldg(j,t)%v)
                                        oldg(j,t)%n=0
                                enddo
                                do j=1,je0
                                        if (oldg(j,t)%n>1) then
                                                do l1=0,3
                                                        m=count(db(oldg(j,t)%v)==l1)
                                                        if (m>0) then
                                                                je1=je1+1
                                                                if (je1>lib_nmax) je1=1
                                                                allocate(oldg(je1,t)%v(m))
                                                                oldg(je1,t)%v=pack(oldg(j,t)%v,db(oldg(j,t)%v)==l1)
                                                                oldg(je1,t)%n=m
                                                        endif
                                                enddo
                                        else
                                                nsg=nsg+1
                                                singletons(nsg)=oldg(j,t)%v(1)
                                        endif
                                        deallocate(oldg(j,t)%v)
                                        oldg(j,t)%n=0
                                enddo
                        endif
                        if (je0<lib_nmax) then
                                js=je0+1
                        else        ! this case is very unlikely
                                js=1
                        endif
                        je0=je1
                enddo
        endif
!        print *, 'tour #', i, ' #(>1, 1, 0)slots=', count(oldg(:,t)%n>1), count(oldg(:,t)%n==1), count(oldg(:,t)%n==0)
enddo
if (js<=je0) then ! last scan to collect the last singletons
        do j=js,je0
                if (oldg(j,t)%n==1) then
                        nsg=nsg+1
                        singletons(nsg)=oldg(j,t)%v(1)
                        deallocate(oldg(j,t)%v)
                        oldg(j,t)%n=0
                endif
        enddo
else
        do j=js,lib_nmax
                if (oldg(j,t)%n==1) then
                        nsg=nsg+1
                        singletons(nsg)=oldg(j,t)%v(1)
                        deallocate(oldg(j,t)%v)
                        oldg(j,t)%n=0
                endif
        enddo
        do j=1, je0
                if (oldg(j,t)%n==1) then
                        nsg=nsg+1
                        singletons(nsg)=oldg(j,t)%v(1)
                        deallocate(oldg(j,t)%v)
                        oldg(j,t)%n=0
                endif
        enddo
endif
print *, '#rounds=', r        
print *, 'final je=', je1
print *, 'oldgs_0=', count(oldg(:,t)%n==0), 'oldgs_1=', count(oldg(:,t)%n==1), 'oldgs_gt1=', count(oldg(:,t)%n>1)
nureadsg=nsg
nureads=count(oldg(:,t)%n>1)
print *, 'nbgr=', nureads
print *, 'nsg=', nureadsg
maxsize=maxval(oldg(:,t)%n)
allocate(ureads(nureads,lhash), depth(nureads,nind), tdepth(nureads), &
        ureadsg(nureadsg,lhash), sgref(nureadsg))
do i=1, nsg
        ureadsg(i,:)=libseqir1(singletons(i),:)
        sgref(i)=xref(singletons(i),1) ! vector sgref shows to whom belong each singleton
enddo
print *, 'counts=', count(sgref==1), count(sgref==2), count(sgref==3)
print *, '#ureads_sg=', nureadsg

depth=0
k=0
do j=1, lib_nmax
        if (oldg(j,t)%n>=2) then
                k=k+1
                l=oldg(j,t)%v(1)
                ureads(k,:)=libseqir1(l,:)
                depth(k,xref(oldg(j,t)%v,1))=xref(oldg(j,t)%v,2)
        endif
        if (allocated(oldg(j,t)%v)) deallocate(oldg(j,t)%v)
enddo
call update_tdepth(nureads,nind)
deallocate(libseqir1, db, oldg, singletons)
print *, ' max depth per read = ', maxval(tdepth), maxloc(tdepth,1)
print *, '         (indiv #1) = ', maxval(depth(:,1)), maxloc(depth(:,1),1)
print *, '         (indiv #2) = ', maxval(depth(:,2)), maxloc(depth(:,2),1)
print *, ''
end subroutine sort_all_unique_raws

subroutine quality_check
implicit none
integer:: props(2)
real:: pA, pB

! STEP 3.1 : Reads QC on over_coverage

print *, '   QUALITY CONTROLS   '
print *, ''
print *, '     (1) Reads QC on over-coverage'
print *, '                          number of slots = ', nureads*nind
nel=count(depth(1:nureads,:)>0)
print *, '                     number of used slots = ', nel
sel=sum(depth)
avg=real(sel)/real(nel)
print *, '                            mean coverage = ', avg
std=0.0
do i=1, nureads
        j=count(depth(i,:)>0)
        allocate(rw2(j))
        rw2=real(pack(depth(i,:),depth(i,:)>0))
        std=std+dot_product(rw2,rw2)
        deallocate(rw2)
enddo
std=std/real(nel-1)
std=sqrt(std)
print *, '              coverage standard deviation = ', std
bounds(1)=avg+3*std
j=nint(bounds(1))
print *, '     coverage upper-bound [real, rounded] = ', bounds(1), j
k=0
l=0
do i=1, nureads
        if (count(depth(i,:)>j)>0) then
                l=l+1
        else
                k=k+1
                depth(k,:)=depth(i,:)
                tdepth(k)=tdepth(i)
                ureads(k,:)=ureads(i,:)
        endif
enddo
! depth(k+1:nureads,:)=0
! tdepth(k+1:nureads)=0
nureads=k
allocate(nindcov(nureads))
do i=1, nureads
        nindcov(i)=count(depth(i,:)>0)
enddo
print *, ' # of reads kept/discarded on coverage QC = ', nureads, l
print *, '          # of reads covering all samples = ', count(nindcov==nind)

! STEP 3.2 : Reads QC on GCcontent

print *, ''
print *, '     (2) Reads QC on GCcontent(%)'
print *, ' Distribution of GC% (singletons reads) :'
allocate(GCout(nureadsg))
call get_GCcontent_stats(ureadsg, nureadsg, lhash, lg, GCout, 3)
m=count(GCout==1)
k=0
do i=1, nureadsg
        if (GCout(i)==1) cycle
        k=k+1
        ureadsg(k,:)=ureadsg(i,:)
        sgref(k)=sgref(i)
enddo
deallocate(GCout)
nureadsg=k
print *, ' '
print *, ' Distribution of GC% (non-singletons reads) :'
allocate(GCout(nureads))
call get_GCcontent_stats(ureads, nureads, lhash, lg, GCout, 3)
n=count(GCout==1)
k=0
do i=1, nureads 
        if (GCout(i)==1) cycle
        k=k+1
        tdepth(k)=tdepth(i)
        depth(k,:)=depth(i,:)
        ureads(k,:)=ureads(i,:)
enddo
nureads=k
deallocate(GCout)
print *, ' Ratio of GC% frequencies [ f(GC%-sg)/f(GC%-non-sg) ] = ', &
        (real(nureads+n)/real(nureadsg+m))*(real(m)/real(n))
print *, '           # of singleton reads kept/discarded on GC% = ', nureadsg, m
print *, '           # of non-singl reads kept/discarded on GC% = ', nureads, n


! STEP 3.3 : Samples QC on #ureads 
print *, ''
print *, '     (3) Sample filtering on number of ureads (including singletons): '
print *, '                   threshold set to = ', minureads
do i=1, nind
        indiv(i)%nur=count(depth(1:nureads,i)>0)+ &
                count(sgref(1:nureadsg)==i)
enddo
print *, '             average #ureads/sample = ', getAvg_real(real(indiv(:)%nur),nind)
m=count(indiv(:)%nur<minureads)
print *, '        number of discarded samples = ', m
if (m==nind) then
        print *, ' ERROR :: no samples left for analysis ! '
        stop
endif
allocate(samples(nind-m))
k=0
do i=1, nind
        if (indiv(i)%nur<minureads) then
                indiv(i)%kept=.false.
                indiv(i)%i=0
        else
                k=k+1
!                tdepth(1:nureads)=tdepth(1:nureads)-depth(1:nureads,i)
                depth(1:nureads,k)=depth(1:nureads,i)
                samples(k)=i
                indiv(i)%i=k
        endif
enddo
nind=k
call update_tdepth(nureads,nind)
print *, ' number of samples left for analysis = ', nind

! STEP 3.4 : Samples QC on check for inter-species
if (interspec_check) then
        print *, ''
        print *, '     (4) Sample inter-specific check : '
        allocate(cpm(nind,nind))
        call compute_cross_proportions_matrix
        k=0;
        do i=1, nind
                props=0
                do j=1, nind
                        if (i==j) cycle
                        pA=cpm(i,j)
                        pB=cpm(j,i)
                        if (pA>.5.and.pB>.5) then
                                props(2)=props(2)+1
                        else
                                props(1)=props(1)+1
                        endif
                enddo
                j=samples(i)
                if (props(1)>props(2)) then
                        print *, '       sample ', indiv(j)%id(:len_trim(indiv(j)%id)), &
                                '  discarded :: ', props
                        indiv(j)%i=0
                        indiv(j)%kept=.false.
                else
                        print *, '          ', indiv(j)%id(:len_trim(indiv(j)%id)), &
                                props
                        k=k+1
        !                tdepth(1:nureads)=tdepth(1:nureads)-depth(1:nureads,i)
                        depth(1:nureads,k)=depth(1:nureads,i)
                        samples(k)=samples(i)
                        indiv(j)%i=k
                endif
        enddo
        call update_tdepth(nureads,nind)
endif

! STEP 3.4 : UPDATE OF UREADS, DEPTH, TDEPTH AND SAMPLES
naddureadsg=0
do i=1, nureads
        if (tdepth(i)==1) naddureadsg=naddureadsg+1
enddo
print *, '   (additional singletons reads after QC = ', naddureadsg, ')'
allocate(addureadsg(naddureadsg,lhash), addsgref(naddureadsg), &
        w2(nind))
! update samples
w2=samples(1:nind)
deallocate(samples)
allocate(samples(nind))
samples=w2
! update ureads (discard unobserved reads) and record new singletons in
! addureadsg, update depth
w2=(/(i,i=1,nind)/)
k=0
j=0
naddureadsg=0
print *, 'nureads=', nureads
do i=1, nureads
        if (tdepth(i)==0) then
                j=j+1
        elseif (tdepth(i)==1) then
                naddureadsg=naddureadsg+1
                addureadsg(naddureadsg,:)=ureads(i,:)
                allocate(w4(1))
                w4=pack(w2,depth(i,1:nind)==1)
                addsgref(naddureadsg)=w4(1)
                deallocate(w4)
        elseif (tdepth(i)>1) then
                k=k+1
                ureads(k,:)=ureads(i,:)
                depth(k,:)=depth(i,:)
        endif
enddo
deallocate(w2)
print *, '   (additional non-singl ureads removed after QC = ', j, ')'
print *, 'k, naddureadsg=', k, naddureadsg
nureads=k
call update_tdepth(nureads,nind)

end subroutine quality_check

subroutine setup_nfemvec
implicit none
integer:: i, j
integer*1:: i1
nfemvec=0
do i1=-128, 127
        i=i1+129
        do j=0, 7
                if (btest(i1,j)) nfemvec(i)=nfemvec(i)+1
        enddo
enddo
end subroutine setup_nfemvec

subroutine setup_GCvec
implicit none
integer:: i, j, k
integer*1:: i1, w(4)
w=(/0,1,0,1/)
GCvec=0
do i1=-128, 127
        i=i1+129
        do j=6, 0, -2
                k=ibits(i1,j,2)
                GCvec(i)=GCvec(i)+w(k+1)
        enddo
enddo
end subroutine setup_GCvec

function get_nfem(v,l) result(n)
implicit none
integer:: l, n, i, j
integer*1:: v(l)
n=0
do i=1, l
        j=v(i)+129
        n=n+nfemvec(j)
enddo
end function get_nfem

function get_GCcontent(v,l,w,lg) result(n)
implicit none
integer:: i, l, lg, w(l)
integer*1:: v(l)
real:: n
w=v+129
w=GCvec(w)
n=real(sum(w))/real(lg)
end function get_GCcontent

function getAvg_real(v,l) result(a)
implicit none
integer:: l
real:: a, v(l)
a=real(sum(v))/real(l)
end function getAvg_real

function getVar_real(v,l) result(s2)
implicit none
integer:: l
real:: a, s2, v(l)
real, allocatable:: w(:)
a=getAvg_real(v,l)
allocate(w(l))
w=real(v)-a
s2=dot_product(w,w)/real(l-1)
deallocate(w)
end function getVar_real

function getStd_real(v,l) result(s)
implicit none
integer:: l
real:: s2, s, v(l)
s2=getVar_real(v,l)
s=s2**(.5)
end function getStd_real

function getNormalBounds_real(v,l,nstd) result(b)
implicit none
integer:: l, nstd
real:: a, s, b(2), v(l)
a=getAvg_real(v,l)
s=getStd_real(v,l)
b(1)=a-nstd*s
b(2)=a+nstd*s
end function getNormalBounds_real

subroutine get_GCcontent_stats(lib, n, l, l0, isout, nstd)
implicit none
integer:: i, j, k, l, l0, n, nstd
integer, allocatable:: w6(:)
integer*1:: lib(n,l), isout(n)
real:: bounds(2)
real, allocatable:: rw2(:)
allocate(w6(l), rw2(n))
do i=1, n
        rw2(i)=get_GCcontent(lib(i,:),l,w6,l0)
enddo
bounds=getNormalBounds_real(rw2,n,nstd)
print *, '            min=', minval(rw2)
print *, '    lower-bound=', bounds(1)
print *, '            avg=', getAvg_real(rw2,n)
print *, '    upper-bound=', bounds(2)
print *, '            max=', maxval(rw2)
print *, '         st-dev=', getStd_real(rw2,n)
print *, '      #outliers=', count(rw2<bounds(1).or.rw2>bounds(2))
print *, '  #not-outliers=', count(rw2>=bounds(1).and.rw2<=bounds(2))
isout=1
forall(i=1:n, rw2(i)>=bounds(1).and.rw2(i)<=bounds(2)) isout(i)=0
end subroutine get_GCcontent_stats


subroutine save_us_files
implicit none
integer:: u1, k1
integer, allocatable:: sgnb(:)
! export data to file for future uses
! first record: #ureads in library | #ureads(inds) | length(hash) | length(seq)
! data (ind) = ureads (ind) | depth (ind) | issg (ind)
allocate(sgnb(nind0))
open(newunit=u1, file='lane_parsing.log')
write(u1,'(a)') '! PLEASE DO NOT EDIT THIS FILE !'
write(u1,*) nureads, lg, lg1, lhash
write(u1,*) count(indiv(:)%kept), nind0
! Singletons
sgnb=0
do i=1, nind0
        if (.not.indiv(i)%kept) cycle
        j=count(sgref==i)+count(addsgref==i)
        indiv(i)%nsg=j
        if (j>0) then
                sgnb(i)=j
                open(newunit=indiv(i)%u,file= &
                        indiv(i)%id(:len_trim(indiv(i)%id))//'.usg',form='binary')
                write(indiv(i)%u) nureadsg+naddureadsg, j, lhash, lg
        endif
enddo
do i=1, nureadsg
        if (.not.indiv(sgref(i))%kept) cycle
        write(indiv(sgref(i))%u) ureadsg(i,:)
enddo
do i=1, naddureadsg
        if (.not.indiv(addsgref(i))%kept) cycle
        write(indiv(addsgref(i))%u) addureadsg(i,:)
enddo
do i=1, nind0
        if (.not.indiv(i)%kept) cycle
        if (indiv(i)%nsg>0) close(indiv(i)%u)
enddo
! Non-singletons
allocate(w6(nureads))
w6=(/(i,i=1,nureads)/)
open(newunit=v, file='reads_counts_per_sample')
k1=0
do i=1, nind0
        if (.not.indiv(i)%kept) then
                write(u1,*) indiv(i)%id(:len_trim(indiv(i)%id)), 0, 0
                cycle
        endif
        k1=k1+1
        h=sum(depth(1:nureads,k1))
        j=count(depth(1:nureads,k1)>0)
        l=count(sgref==i)
        write(u1,*) indiv(i)%id(:len_trim(indiv(i)%id)), 1, j
        allocate(ureads1(j,lhash), indiv(i)%libp(j), indiv(i)%libd(j))
        do k=1, lhash
                ureads1(:,k)=pack(ureads(1:nureads,k),depth(1:nureads,k1)>0)
        enddo
        indiv(i)%libd=pack(depth(1:nureads,k1),depth(1:nureads,k1)>0)
        indiv(i)%libp=pack(w6,depth(1:nureads,k1)>0)
        open(newunit=u,file=indiv(i)%id(:len_trim(indiv(i)%id))//'.us',form='binary')
        write(u) nureads, j, lhash, l, lg
        do k=1,j
                write(u) ureads1(k,:), indiv(i)%libd(k), indiv(i)%libp(k)
        enddo
        close(u)
        write(v,'(a,5i12,i5)') indiv(i)%id(:len_trim(indiv(i)%id)), h+l, j+l, h, j, l, lg
        deallocate(ureads1)
enddo
close(v)
deallocate(sgref,w6)
write(u1,'(a)') '! PLEASE DO NOT EDIT THIS FILE !'
close(u1)
end subroutine save_us_files

subroutine print_section_separation
print *, ''
print '(a)', '----------------------------------------------------------------------------------------------'
end subroutine print_section_separation

subroutine print_section_title(title)
implicit none
integer:: i, j, k(2)
character*50:: title, dash, pad1
dash='----------------------'
pad1='                                                  '
title=trim(adjustl(title))
i=len_trim(title)
j=50-i
if (mod(j,2)==0) then
        k(1)=j/2
        k(2)=k(1)
else
        k(1)=(j-1)/2
        k(2)=k(1)+1
endif
call print_section_separation
print '(a,a,a,a,a)', dash(:len_trim(dash)), pad1(:k(1)), title(:len_trim(title)), &
         pad1(:k(2)), dash(:len_trim(dash))
print *, ''
end subroutine print_section_title

subroutine compute_cross_proportions_matrix
implicit none
integer:: i, j, naub, na
real:: cp(2)
cpm=0.0
do i=1, nind
        do j=1, nind
                if (i==j) cycle
                naub=sum(depth(1:nureads,i),MASK=depth(1:nureads,j)>0)
                na=sum(depth(1:nureads,i))
                cpm(i,j)=real(naub)/real(na)
        enddo
enddo
end subroutine compute_cross_proportions_matrix

subroutine load_binary_files
implicit none
integer:: u1, iskept, nb, nind0p
integer*1, allocatable:: xrf(:)
logical:: fex
inquire(file='lane_parsing.log', exist=fex)
if (.not.fex) then
        print *, ' ERROR :: No file lane_parsing.log found in the current directory'
        stop
endif
open(newunit=u1, file='lane_parsing.log')
read(u1,*) ! skip header line
read(u1,*) nureads, lg, lg1, lhash
read(u1,*) nind, nind0p
print *, '    # of individuals parsed in lanes = ', nind0p
print *, '  # of individuals kept for analyses = ', nind
print *, '                    good read length = ', lg
print *, '                  hashed read length = ', lhash
print *, '               un-hashed read length = ', lg1
print *, '              # of ureads in library = ', nureads
allocate(samples(nind))
k=0
do i=1, nind0p
        read(u1,*) junk20, iskept, nb
        if (iskept==0) then
                indiv(i)%i=0
                indiv(i)%kept=.false.
        else
                k=k+1
                indiv(i)%i=k
                samples(k)=i
                indiv(i)%kept=.true.
                indiv(i)%nur=nb
                open(newunit=indiv(i)%u,file=indiv(i)%id(:len_trim(indiv(i)%id))//'.us',form='binary')
                read(indiv(i)%u) nureads, nb, lhash, nsg, lg
                if (trim(adjustl(junk20))/=indiv(i)%id) then
                        print *, ' WARNING : no match between IDs :', &
                                trim(adjustl(junk20)), indiv(i)%id
                endif
        endif
enddo
nureads0=sum(indiv(:)%nur)
print *, '            total #ureads in library = ', nureads
print *, '              total #ureads in files = ', nureads0
allocate(ureads(nureads,lhash), w2i1(lhash), depth(nureads,nind), &
        tdepth(nureads), xrf(nureads))
xrf=0
k=0
depth=0
tdepth=0
do i=1, nind0
        if (.not.indiv(i)%kept) cycle
        k=k+1
        print *,  '       .... loading sample ', &
                indiv(i)%id(:len_trim(indiv(i)%id)), ' ....'
        allocate(indiv(i)%libp(indiv(i)%nur), indiv(i)%libd(indiv(i)%nur))
        do j=1, indiv(i)%nur
                read(indiv(i)%u) w2i1, indiv(i)%libd(j), l
                if (xrf(l)==0) then
                        ureads(l,:)=w2i1
                        xrf(l)=1
                endif
                indiv(i)%libp(j)=l
        enddo
        depth(indiv(i)%libp,k)=indiv(i)%libd
!        tdepth(indiv(i)%libp)=tdepth(indiv(i)%libp)+indiv(i)%libd
enddo   
deallocate(w2i1, xrf)
call update_tdepth(nureads,nind)
end subroutine load_binary_files

subroutine apply_depth_filters
implicit none
integer:: k, i
k=0
do i=1, nureads
        if (tdepth(i)<tdepthmin) cycle
        if (maxval(depth(i,:))<depthmin) cycle
        k=k+1
        tdepth(k)=tdepth(i)
        depth(k,:)=depth(i,:)
        ureads(k,:)=ureads(i,:)
enddo
nureads=k
print *, '#ureads (>= tdepthmin, >= depthmin)=', nureads
end subroutine apply_depth_filters

subroutine sex_identification
implicit none
integer:: lhashpat
character*12:: num
character*10:: header
call set_nfemvec
call cpu_time(t0)
lhashpat=nind
print *, '              lenght of hash of patterns = ', lhashpat
npats=nureads
allocate(patterns(npats,lhashpat), w4(npats))
w4=(/(i,i=1,nureads)/)
patterns=0
k=0
do i=1, lhashpat
                forall(l=1:npats, depth(w4(l),i)>0) patterns(l,i)=1
enddo
deallocate(w4)
t=1
js=1
je0=1
r=0
nsg=0
allocate(oldg(lib_nmax,1), singletons(sg_nmax))
allocate(oldg(1,t)%v(npats))
oldg(:,t)%n=0
oldg(1,t)%v=(/(i,i=1,npats)/)
oldg(1,t)%n=npats 
do i=1,lhashpat
        je1=je0
        if (js<=je0) then
                do j=js,je0
                        if (oldg(j,t)%n>1) then
                                do l1=0, 1
                                        m=count(patterns(oldg(j,t)%v,i)==l1)
                                        if (m>0) then
                                                je1=je1+1
                                                if (je1>lib_nmax) je1=1
                                                allocate(oldg(je1,t)%v(m))
                                                oldg(je1,t)%v=pack(oldg(j,t)%v,patterns(oldg(j,t)%v,i)==l1)
                                                oldg(je1,t)%n=m
                                        endif
                                enddo
                        else
                                nsg=nsg+1
                                singletons(nsg)=oldg(j,t)%v(1)
                        endif
                        deallocate(oldg(j,t)%v)
                        oldg(j,t)%n=0
                enddo
        else
                r=r+1
                do j=js,lib_nmax
                        if (oldg(j,t)%n>1) then
                                do l1=0, 1
                                        m=count(patterns(oldg(j,t)%v,i)==l1)
                                        if (m>0) then
                                                je1=je1+1
                                                if (je1>lib_nmax) je1=1
                                                allocate(oldg(je1,t)%v(m))
                                                oldg(je1,t)%v=pack(oldg(j,t)%v,patterns(oldg(j,t)%v,i)==l1)
                                                oldg(je1,t)%n=m
                                        endif
                                enddo
                        else
                                nsg=nsg+1
                                singletons(nsg)=oldg(j,t)%v(1)
                        endif
                        deallocate(oldg(j,t)%v)
                        oldg(j,t)%n=0
                enddo
                do j=1,je0
                        if (oldg(j,t)%n>1) then
                                do l1=0, 1
                                        m=count(patterns(oldg(j,t)%v,i)==l1)
                                        if (m>0) then
                                                je1=je1+1
                                                if (je1>lib_nmax) je1=1
                                                allocate(oldg(je1,t)%v(m))
                                                oldg(je1,t)%v=pack(oldg(j,t)%v,patterns(oldg(j,t)%v,i)==l1)
                                                oldg(je1,t)%n=m
                                        endif
                                enddo
                        else
                                nsg=nsg+1
                                singletons(nsg)=oldg(j,t)%v(1)
                        endif
                        deallocate(oldg(j,t)%v)
                        oldg(j,t)%n=0
                enddo
        endif
        if (je0<lib_nmax) then
                js=je0+1
        else        ! this case is very unlikely
                js=1
        endif
        je0=je1
enddo
if (js<=je0) then ! last scan to collect the last singletons
        do j=js,je0
                if (oldg(j,t)%n==1) then
                        nsg=nsg+1
                        singletons(nsg)=oldg(j,t)%v(1)
                        deallocate(oldg(j,t)%v)
                        oldg(j,t)%n=0
                endif
        enddo
else
        do j=js,lib_nmax
                if (oldg(j,t)%n==1) then
                        nsg=nsg+1
                        singletons(nsg)=oldg(j,t)%v(1)
                        deallocate(oldg(j,t)%v)
                        oldg(j,t)%n=0
                endif
        enddo
        do j=1, je0
                if (oldg(j,t)%n==1) then
                        nsg=nsg+1
                        singletons(nsg)=oldg(j,t)%v(1)
                        deallocate(oldg(j,t)%v)
                        oldg(j,t)%n=0
                endif
        enddo
endif
print *, '#rounds=', r        
print *, 'final je=', je1
print *, 'oldgs_0=', count(oldg(:,t)%n==0), 'oldgs_1=', count(oldg(:,t)%n==1), 'oldgs_gt1=', count(oldg(:,t)%n>1)
nupats=count(oldg(:,t)%n>1)+nsg
print *, 'nbgr=', nupats
print *, 'nsg=', nsg
maxsize=maxval(oldg(:,t)%n)
print *, 'maxsize=', maxsize
allocate(nbocc(nupats), upatterns(nupats,lhashpat), sexv(lhashpat)) !!! TODO check
!! everything in order that sexv has length nind (lhashpat because patterns are not
!! coded)
nbocc=0
do i=1, nsg
        upatterns(i,:)=patterns(singletons(i),:)
        nbocc(i)=1
enddo
k=nsg
do j=1, lib_nmax
        if (oldg(j,t)%n>=2) then
                k=k+1
                l=oldg(j,t)%v(1)
                upatterns(k,:)=patterns(l,:)
                nbocc(k)=oldg(j,t)%n
        endif
        if (allocated(oldg(j,t)%v)) deallocate(oldg(j,t)%v)
enddo
deallocate(patterns, oldg, singletons)
print *, ' number of singleton patterns =', count(nbocc==1)
print *, ' number of errors =', count(nbocc==0)
call cpu_time(t1)
print *, ' elapsed time for pattern sorting =', t1-t0, ' seconds'

! Find the "most likely" sex distribution pattern, considering those with
! occurrence >= threstopocc
call cpu_time(t0)
call set_log_pasctri
threstopocc=nint(1E-12*real(npats))
print *, ' -- Threshold for top occurences =', threstopocc
ntopocc=count(nbocc>threstopocc)
allocate(w2(nupats), topocc(ntopocc,3), scores(ntopocc,3))
w2=nbocc
maxsc1=0.0
maxsc2=0.0
print *, ' -- Patterns with occurence>.05% (rank, #occ, #hemi, score):: #=', ntopocc
ntopocc=0
do i=1, nupats
        j=maxloc(w2,1)
        k=nbocc(j)
        if (k<=threstopocc) exit
        ntopocc=ntopocc+1
        l=get_nfem(upatterns(j,:),lhashpat)
        topocc(ntopocc,:)=(/k,j,l/)
        scores(ntopocc,2)=2.d0**-lpt(l+1)
        scores(ntopocc,3)=scores(ntopocc,2)*real(nbocc(j))
        if (scores(ntopocc,3)>maxsc1) then
                maxsc1=scores(ntopocc,3)
                best1=j
        endif
        if (scores(ntopocc,3)>maxsc2.and.scores(ntopocc,3)<maxsc1) then
                maxsc2=scores(ntopocc,3)
                best2=j
        endif
        w2(j)=0
enddo
deallocate(w2)
call cpu_time(t1)
print *, ' elapsed time for obtaining sex scores =', t1-t0, ' seconds'
sexv=upatterns(best1,:)
print *, ' -- Top-2 highest-score patterns (score, #occ, #hemi)='
print *, '     ', maxsc1, nbocc(best1), get_nfem(upatterns(best1,:), lhashpat)
print *, '     ', maxsc2, nbocc(best2), get_nfem(upatterns(best2,:), lhashpat)
print *, ' -- Confidence score for prediction = ', 1.d0-maxsc2/maxsc1
if (optim_opt==1) then
        print *, ' -- Confidence score optimization ::'
        j=nint(1E-2*nsexminpc*real(nind))
        k=nind-j
        maxsc1=0.d0
        do i=j, k
                maxsc11=0.d0; maxsc21=0.d0
                sexrate=real(i)/real(nind)
                call set_log_pasctri
                do l=1, ntopocc
                        scores(l,2)=2.d0**-lpt(topocc(l,3)+1)
                        scores(l,3)=scores(l,1)*scores(l,2)*real(topocc(l,1))
                        if (scores(l,3)>maxsc11) then
                                maxsc11=scores(l,3)
                                best11=topocc(l,2)
                        endif
                        if (scores(l,3)>maxsc21.and.scores(l,3)<maxsc11) then
                                maxsc21=scores(l,3)
                        endif
                enddo
                print '(a,2f7.2,2f12.3)', '    ', sexrate, 1.d0-(maxsc21/maxsc11), maxsc11, maxsc21 
                if (maxsc1<1.d0-(maxsc21/maxsc11)) then
                        maxsc1=1.d0-(maxsc21/maxsc11)
                        best1=best11
                        bestsexrate=sexrate
                endif
        enddo                
        sexv=upatterns(best1,:)
        print *, ' -- Optimized sex-rate =', bestsexrate
endif
open(newunit=u, file='sex_assignment')
k=0
do i=1, nind0
        if (.not.(indiv(i)%kept)) cycle
        k=k+1
        write(u,'(a,a,a)') indiv(i)%id(:len_trim(indiv(i)%id)), ' ', sexc(sexv(k)+1)
enddo
close(u)
! spot the sex-linked reads
j=count(sexv==0)
k=count(sexv==1)
allocate(indsM(j), indsF(k), w2(nind), isW(nureads))
w2=(/(i,i=1, nind)/)
indsM=pack(w2,sexv==0)
indsF=pack(w2,sexv==1)
deallocate(w2)
isW=.false.
do i=1, nureads
        if (minval(depth(i,indsF))>0.and.sum(depth(i,indsM))==0) &
                isW(i)=.true.
enddo
deallocate(indsM, indsF)
! if needed, output sex-linked reads
if (save_opt) then
        open(newunit=u, file='likelyW.fasta')
        open(newunit=v, file='likelyW.depth')
        allocate(seq0(lg))
        write(fmt1,'(a,i4,a)') '(',nind,'i8)'
        write(fmt2,'(a,i4,a)') '(', lg, 'a1)'
        write(header,'(a10)') '>sequence_'
        do i=1, nureads
                if (isW(i)) then
                        write(num,'(i12)') i
                        write(u,'(a)') header//trim(adjustl(num))
                        seq0=seqir2reads(ureads(i,:),lhash,lg)
                        write(u,fmt2) seq0
                        write(v,fmt1) depth(i,1:nind)
                endif
        enddo
        deallocate(seq0)
        close(u)
        close(v)
endif
! update ureads, depth, tdepth by removing sex-linked reads
k=0
do i=1, nureads
        if (isW(i)) cycle
        k=k+1
        ureads(k,:)=ureads(i,:)
        depth(k,:)=depth(i,:)
enddo
nureads=count(.not.isW)
nlib=nureads
call update_tdepth(nureads, nind)
print *, '#ureads in library:', nureads
end subroutine sex_identification

subroutine set_nfemvec
implicit none
integer:: i, j
integer*1:: i1
nfemvec=0
do i1=-128, 127
        i=i1+129
        do j=0, 7
                if (btest(i1,j)) nfemvec(i)=nfemvec(i)+1
        enddo
enddo
end subroutine set_nfemvec

subroutine setup_translation_vec
implicit none
integer:: i, n, j, k, l
integer*1:: i1, alphanum(4)
character*1:: alphabet(4)
alphabet=(/'A','C','T','G'/)
alphanum=(/0,1,2,3/)
n=256
allocate(seqirvec(n,4), seqivec(n,4))
do i1=-128,127
        i=i1+129
        l=0
        do j=6,0,-2
                l=l+1
                k=ibits(i1,j,2)+1
                seqirvec(i,l)=alphabet(k)
                seqivec(i,l)=alphanum(k)
        enddo
enddo
end subroutine setup_translation_vec

function seqir2reads(s,lh,l) result(r)
implicit none
integer:: lh, l, i, j, k
integer*1:: s(lh)
character*1:: r(l), last(4)
do i=1, lh-1
        j=s(i)+129
        k=(i-1)*4
        r(k+1:k+4)=seqirvec(j,:)
enddo
last=seqirvec(s(lh)+129,:)
do i=k+5,l
        r(i)=last(i-k-4)
enddo
end function seqir2reads

subroutine set_log_pasctri
implicit none
integer:: i, j
real*8:: r, two, nr, jr
real*8, allocatable:: v(:), w(:)
two=2.d0
allocate(w(1),v(1))
if (not(allocated(lpt))) allocate(lpt(nind+1))
w(1)=0.d0
do j=1, nind
        deallocate(v)
        allocate(v(j+1))
        v(1)=real(j)
        v(j+1)=real(j)
        do i=2, j
                r=real(j)/real(1+j-i)
                v(i)=1+w(i)-log(r)/log(two)
        enddo
        deallocate(w)
        allocate(w(j+1))
        w=v
enddo
if (sexrate==0.5) then
        lpt=w
else
        nr=real(nind,8)
        do j=0, nind
                r=w(j+1)
                jr=real(j,8)
                lpt(j+1)=r-nr-jr*log(sexrate)/log(two)-(nr-jr)*log(1.d0-sexrate)/log(two)
        enddo
endif
end subroutine set_log_pasctri

subroutine update_tdepth(n,p)
implicit none
integer:: n, p, i
do i=1, n
        tdepth(i)=sum(depth(i,1:p))
enddo
end subroutine update_tdepth

subroutine save_library
implicit none
character*12:: num
character*10:: header
open(newunit=u, file='library.fasta')
open(newunit=v, file='library.depth')
write(fmt1,'(a,i3,a)') '(',nind,'i8)'
write(fmt2,'(a,i3,a)') '(', lg, 'a1)'
write(header,'(a10)') '>sequence_'
do i=1, nureads
        write(num,'(i12)') i
        write(u,'(a)') header//trim(adjustl(num))
        seq0=seqir2reads(ureads(i,:),lhash,lg)
        write(u,fmt2) seq0
        write(v,fmt1) depth(i,:)
enddo
close(u)
close(v)
end subroutine save_library

subroutine compar_reads
implicit none
! TODO update the group_countdown function in order to use the grouping function
!       the whash2, khash2, kkhash2 and lhash2 variables refer to that former
!       group_countdown function
integer:: h, i, j, k, l, m, n
integer, allocatable:: w1(:), w2(:), p(:), w5(:,:)
integer*8:: nbc
real:: sav

print *, '# of allowed differences =', nbdif
print *, 'maximum size of indels (pb) =', nbpb
allocate(tags(ntagmax,3))
ntag=0; ntag0=0; nbgr=0
do ndi= 1, nbdif
        call splitting(ndi)
        do sp=1, nsplit
                l=split(sp,1)%n
                lhash2=ceiling(real(l)/real(khash2))
                allocate(libseqir(nlib,lhash2), w5(nlib,l))
                w5=libseqi(:,split(sp,1)%v)
                libseqir=hash_library(w5,nlib,l,lhash2)
                call initialize_newg(nbgr)
                call group_countdown(libseqir,nlib,lhash2,nbgr,maxsize,nbc,sav)
                deallocate(w5)
                l=split(sp,2)%n
                allocate(w5(nlib,l))
                w5=libseqi(:,split(sp,2)%v)
                call find_snps_in_groups(nbgr,w5,nlib,l,ndi)
                print '(a,i2,a,i1,a,i8,i12,i9,f9.4,i8)', '       split #', sp, ' : (', ndi, ') : ', &
                        maxsize, nbc, nbgr, 100.d0*sav, ntag-ntag0
                deallocate(w5, libseqir)
                ntag0=ntag
        enddo
enddo
print '(a, i12)', 'final # of tags= ', ntag
do i=1, nbdif
        print '(a,i1,a,i8)', '      #tags with ', i,' diff =', count(tags(1:ntag,3)==i)
enddo
! resizing of tags
allocate(w3(ntag,3))
w3=tags(1:ntag,:)
call move_alloc(w3,tags)

if (nbpb>0) then
        ! launch indels search
        allocate(indels(4500000,4), liborigin(nlib*2,2))
        liborigin(1:nlib,1)=0; liborigin(nlib+1:nlib*2,1)=1
        liborigin(1:nlib,2)=(/(i,i=1,nlib)/)
        liborigin(nlib+1:nlib*2,2)=liborigin(1:nlib,2)
        nindel=0; nindel0=0; j=nint(real(lg1)/2.d0)
        do npi=1, nbpb
                !forward search
                k=j-npi; l=k
                lhash2=ceiling(real(l)/real(khash2))
                allocate(libseqir(nlib,lhash2), w5(nlib,l))
                w5=libseqi(:,1:l)
                libseqir=hash_library(w5,nlib,l,lhash2)
                call initialize_newg(nbgr)
                call group_countdown(libseqir,nlib,lhash2,nbgr,maxsize,nbc,sav)
                deallocate(w5)
                l=lg1-k
                allocate(w5(nlib,l),w1(l))
                w5=libseqi(:,k+1:lg1)
                w1=(/(i, i=k+1,lg1)/)
                call find_indels_forward(nbgr,w5,nlib,l,npi,w1)
                print '(a,i1,a,i8,i12,i9,f9.4,i8)', '  forward search for ', npi,'-indels : ', & 
                        maxsize, nbc, nbgr, 100.d0*sav, nindel-nindel0
                deallocate(libseqir,w5,w1)
                nindel0=nindel
                !backward search
                l=lg1-j-npi
                lhash2=ceiling(real(l)/real(khash2))
                allocate(libseqir(nlib*2,lhash2), w5(nlib*2,l))
                w5(1:nlib,:)=libseqi(:,j+1:lg1-npi)
                w5(nlib+1:nlib*2,:)=libseqi(:,j+1+npi:lg1)
                libseqir=hash_library(w5,nlib*2,l,lhash2)
                call initialize_newg(nbgr)
                call group_countdown(libseqir,nlib*2,lhash2,nbgr,maxsize,nbc,sav)
                deallocate(w5)
                l=j+npi
                allocate(w5(nlib*2,l),w1(l))
                w5(1:nlib,:)=libseqi(:,1:j+npi)
                w5(nlib+1:nlib*2,:)=libseqi(:,1:j+npi)
                w1=(/(i, i=1,j+npi)/)
                call find_indels_backward(nbgr,w5,nlib*2,l,npi,w1,j)
                print '(a,i1,a,i8,i12,i9,f9.4,i8)', ' backward search for ', npi,'-indels : ', & 
                        maxsize, nbc, nbgr, 100.d0*sav, nindel-nindel0
                deallocate(libseqir,w5,w1)
                nindel0=nindel
        enddo
        ! resizing of indels
        allocate(w3(nindel,3))
        w3=indels(1:nindel,:)
        call move_alloc(w3,indels)
endif
do i=1, nbpb
        print '(a,i1,a,i8)', '#indel-tags with ', i,' pb =', count(indels(1:nindel,3)==i)
enddo
print '(a,i8)', ' Total number of indel-tags = ', nindel
if (save_opt) then
        open(newunit=j, file='tags_library')
        do i=1, ntag
                write(j,'(2i7,i2)') tags(i,1), tags(i,2), tags(i,3)
        enddo
        close(j)
        open(newunit=j, file='indels_library')
        do i=1, nindel
                write(j,'(2i7,i2,i3)') indels(i,:)
        enddo
        close(j)
endif
end subroutine compar_reads

subroutine splitting(nd)
implicit none
integer:: i, j, k, l, nel, c, c1, c2, nd
integer, allocatable:: table(:,:), w1(:,:)
if (allocated(split)) then
        deallocate(split)
endif
nsplit=2*(2**nd-1)
allocate(table(nsplit,nd+1), split(nsplit,2), p(nd), w1(nd+1,3))
table=0 ! filling only until nbdif = 3 at maximum
do i=1, nd+1
        table(i,i)=1
enddo
if (nd>1) then
        do i=1, nd+1
                table(nd+1+i,:)=abs(1-table(i,:))
        enddo
        c=8
        if (nd>2) then
                do i=2,nd+1
                        do j=1,i-1
                                c=c+1
                                table(c,i)=1
                                table(c,j)=1
                        enddo
                enddo
        endif
endif
do i=1, nd
        p(i)=nint(real(i)*real(lg1)/real(nd+1))
enddo
w1(1,1)=1
w1(1:nd,2)=p
w1(2:nd+1,1)=p+1
w1(nd+1,2)=lg1
w1(:,3)=w1(:,2)-w1(:,1)+1
do i=1, nsplit
        nel=0
        do j=1, nd+1
                if (table(i,j)==0) cycle
                nel=nel+w1(j,3)
        enddo
        allocate(split(i,1)%v(nel), split(i,2)%v(lg1-nel))
        split(i,1)%n=nel; split(i,2)%n=lg1-nel
        c1=0; c2=0
        do j=1, nd+1
                if (table(i,j)==0) then
                        split(i,2)%v(c1+1:c1+w1(j,3))=(/(k, k=w1(j,1),w1(j,2))/)
                        c1=c1+w1(j,3)
                else
                        split(i,1)%v(c2+1:c2+w1(j,3))=(/(k, k=w1(j,1),w1(j,2))/)
                        c2=c2+w1(j,3)
                endif
        enddo
enddo
deallocate(table,p,w1)
end subroutine splitting

subroutine find_indels_backward(ng,lib,n,m,nd,pos,st)
implicit none
integer:: ng, n, m, lib(n,m), h, i, j, k ,l, nd, pos(m), sd, f, st, A, B
integer, allocatable:: w1(:), w2(:,:)
! We start in position st as the comparisons are supposed to be already done on further
! position (i.e. they were grouped together)
! The insert carrier is denoted B and the deletion carrier is denoted A
do i=1, ng
        allocate(w1(newg(i)%n))
        w1=newg(i)%v
        do j=2, newg(i)%n
                do k=1, j-1
                        if (liborigin(w1(j),1)==liborigin(w1(k),1)) cycle
                        if (liborigin(w1(j),2)==liborigin(w1(k),2)) cycle
                        if (liborigin(w1(j),1)==0) then
                                A=w1(j)
                                B=w1(k)
                        else
                                A=w1(k)
                                B=w1(j)
                        endif
                        do h=st, 1, -1
                                if (lib(A,h)/=lib(B,h+nd)) exit
                        enddo
                        if (h==0) cycle
                        sd=count(lib(A,1:h)/=lib(B,1:h))
                        if (sd==0) then
                                f=0
                                do l=2, h
                                        if (lib(A,l)/=lib(A,l-1)) then
                                                f=1
                                                exit
                                        endif
                                enddo
                                if (f==1) then
                                        ! check if that indel was not yet recorded (the smallest = the best)
                                        l=count(indels(1:nindel,1)==liborigin(A,2).and.indels(1:nindel,2)==liborigin(B,2))
                                        if (l==0) then
                                                nindel=nindel+1
                                                indels(nindel,1:4)=(/liborigin(A,2),liborigin(B,2),nd,pos(h+1)/)
                                        endif
                                endif
                        endif
                enddo
        enddo
        deallocate(w1)
enddo
end subroutine find_indels_backward

subroutine find_indels_forward(ng,lib,n,m,nd,pos)
implicit none
integer:: ng, n, m, lib(n,m), i, j, k ,l, c, h, nd, pos(m), sd, f
integer, allocatable:: w1(:), w2(:,:)
do i=1, ng
        allocate(w1(newg(i)%n))
        w1=newg(i)%v
        do j=2, newg(i)%n
                do k=1, j-1
                        do h=1, m-nd
                                if (lib(w1(j),h)/=lib(w1(k),h)) exit
                        enddo
                        if (h>=m-nd) cycle
                        c=h
                        sd=count(lib(w1(j),c+nd:m)/=lib(w1(k),c:m-nd))
                        if (sd==0) then
                                !check if not a mononorphic sequence
                                f=0
                                do  l=c+nd+1, m
                                        if (lib(w1(j),l)/=lib(w1(j),l-1)) then
                                                f=1
                                                exit
                                        endif
                                enddo
                                if (f==1) then
                                        !check if that indel exist (rule the smallest the best)
                                        l=count(indels(1:nindel,1)==w1(k).and.indels(1:nindel,2)==w1(j))
                                        if (l==0) then
                                                nindel=nindel+1
                                                indels(nindel,1:4)=(/w1(k),w1(j),nd,pos(c)/)
                                        endif
                                endif
                        endif
                        sd=count(lib(w1(k),c+nd:m)/=lib(w1(j),c:m-nd))
                        if (sd==0) then
                                !check if not a mononorphic sequence
                                f=0
                                do  l=c+nd+1, m
                                        if (lib(w1(k),l)/=lib(w1(k),l-1)) then
                                                f=1
                                                exit
                                        endif
                                enddo
                                if (f==1) then
                                        l=count(indels(1:nindel,1)==w1(j).and.indels(1:nindel,2)==w1(k))
                                        if (l==0) then
                                                nindel=nindel+1
                                                indels(nindel,1:4)=(/w1(j),w1(k),nd,pos(c)/)
                                        endif
                                endif
                        endif
                enddo
        enddo
        deallocate(w1)
enddo
end subroutine find_indels_forward


subroutine find_snps_in_groups(ng,lib,n,m,nd)
implicit none
integer:: ng, n, m, lib(n,m), i, j, k ,l, c, h, nd
integer, allocatable:: w1(:), w2(:,:)
do i=1, ng
        allocate(w1(newg(i)%n))
        w1=newg(i)%v
        do j=2, newg(i)%n
                do k=1, j-1
                        c=0
                        do h=1, m
                                if (lib(w1(j),h)/=lib(w1(k),h)) c=c+1
                                if (c>nd) exit
                        enddo
                        if (c/=nd) cycle
                        if (ntag==ntagmax) then
                                l=ntagmax+5E5
                                allocate(w2(l,3))
                                w2(1:ntag,1:3)=tags(1:ntag,1:3)
                                deallocate(tags)
                                allocate(tags(l,3))
                                tags(1:ntag,1:3)=w2(1:ntag,1:3)
                                deallocate(w2)
                                ntagmax=l
                        endif
                        ntag=ntag+1
                        tags(ntag,1:3)=(/w1(k),w1(j),c/)
                enddo
        enddo
        deallocate(w1)
enddo
end subroutine find_snps_in_groups

function hash_library(libin,n,m,mb) result(lib)
implicit none
integer:: n, m, mb, libin(n,m), lib(n,mb), i, j, k, l
do i=1,mb-1
        k=i*khash2
        j=k-khash2+1
        do l=1, n
                lib(l,i)=1+dot_product(libin(l,j:k),whash2)
        enddo
enddo
j=k+1
k=m-j+1
do l=1, n
        lib(l,mb)=1+dot_product(libin(l,j:m),whash2(1:k))
enddo
end function hash_library

subroutine group_countdown(lib,n,m,ng,maxsize,nbc,sav)
!lib must be a table of values hashed using khash2, whash2, lhash2, kkhash2
implicit none
integer:: n, m, lib(n,m), ng, i, j, k, l, c, maxsize
integer, allocatable:: w1(:), xrf(:)
integer*8:: maxv, jadd1, jadd2, nbc
real:: sav
type(icell), allocatable:: gtmp(:)
maxv=m*4**khash2+1
allocate(xrf(maxv), gtmp(n))
xrf=0; ng=0
do i=1, n
        gtmp(i)%n=0
        jadd1=1
        do j=1,m
                jadd1=jadd1+mod(kkhash2,lib(i,j))
        enddo
        if (xrf(jadd1)==0) then
                ng=ng+1
                xrf(jadd1)=i
                allocate(gtmp(i)%v(1))
                gtmp(i)%v(1)=i
                gtmp(i)%n=1
                gtmp(i)%i=ng
        else
                jadd2=jadd1
                do
                        j=xrf(jadd2)
                        if (j>0) then
                                c=0; l=0
                                do
                                        l=l+1
                                        if (lib(i,l)/=lib(j,l)) then
                                                c=c+1
                                                exit
                                        endif
                                        if (l==m) exit
                                enddo
                                if (c==0) then
                                        l=gtmp(j)%n
                                        allocate(w1(l+1))
                                        w1(1:l)=gtmp(j)%v
                                        w1(l+1)=i
                                        deallocate(gtmp(j)%v)
                                        allocate(gtmp(j)%v(l+1))
                                        gtmp(j)%v=w1
                                        gtmp(j)%n=l+1
                                        deallocate(w1)
                                        exit
                                endif
                        else
                                ng=ng+1
                                xrf(jadd2)=i
                                gtmp(i)%n=1
                                gtmp(i)%i=i
                                allocate(gtmp(i)%v(1))
                                gtmp(i)%v(1)=i
                                exit
                        endif
                        jadd2=jadd2+1
                enddo
        endif
enddo
ng=count(gtmp(:)%n>1)
allocate(newg(ng))
j=0; maxsize=0; nbc=0
do i=1, n
        l=gtmp(i)%n
        if (l<=1) cycle
        if (l>maxsize) maxsize=l
        j=j+1
        allocate(newg(j)%v(l))
        newg(j)%v=gtmp(i)%v
        newg(j)%n=l
        newg(j)%i=gtmp(i)%i
        nbc=nbc+newg(j)%n*(newg(j)%n-1)/2
enddo
do i=1, n
        if (gtmp(i)%n>0) then
                deallocate(gtmp(i)%v) ! proper cleaning of gtmp
                gtmp(i)%n=0
        endif
enddo
deallocate(xrf,gtmp)
sav=1.d0-2.d0*real(nbc)/(real(n)*real(n-1))
end subroutine group_countdown

subroutine initialize_newg(n)
implicit none
integer:: i, n
do i=1, n
        if (allocated(newg(i)%v)) deallocate(newg(i)%v)
enddo
if (allocated(newg)) deallocate(newg)
end subroutine initialize_newg

subroutine get_libseqi
implicit none
integer:: i, j, k, l, n
allocate(libseqi(nlib,lg1))
do i=1, nureads
        l=-3
        do j=1, lhash
                l=l+khash
                n=ureads(i,j)+129
                libseqi(i,l:l+khash-1)=seqivec(n,:)
        enddo
enddo
end subroutine get_libseqi

subroutine biAssemble2
! TODO  Re-write this part in order to enhance computational performances (so
! far, this part is based on the Matlab prototype). When translating to Fortran,
! one (small) bug was found: the tags library may list several times a
! reciprocal pair; in such case, that pair was not considered as unique to the
! two reads in the pair and subsequently not called.
implicit none
integer:: i, j, k, l, ne, ne2, ks, nnsi, nr, maxpairs, v, nm6
nnsi=nbdif+nbpb
allocate(recip(nnsi),redunds(nureads,nnsi))
k=0
do i=1, nbdif
        k=k+1
        call get_reciprocals(i,tags,ntag,3,k)
        print *, '   ',i, ' SNPs, #reciprocals = ', &
                recip(k)%n
enddo
do i=1, nbpb
        k=k+1
        call get_reciprocals(i,indels(:,1:3),nindel,3,k)
        print *, '   ', i, ' indels, #reciprocals = ', &
                recip(k)%n
enddo
redunds=0
do k=1, nnsi
        do l=1, recip(k)%n
                i=recip(k)%v(l,1)
                j=recip(k)%v(l,2)
                redunds(i,k)=l
                redunds(j,k)=l
        enddo
enddo
maxpairs=sum(recip(:)%n)
print *, '       maximum # of pairs = ', maxpairs
allocate(pairs(maxpairs,6))
np=0
do l=1, nureads
        ne=count(redunds(l,:)>0)
        if (ne==0) cycle
        if (ne>1) then
                call gather_web(l, nureads, nnsi)
        elseif (ne==1) then
                k=maxloc(redunds(l,:),1)
                i=recip(k)%v(redunds(l,k),1)
                j=recip(k)%v(redunds(l,k),2)
                if (i==l) then
                        ne2=count(redunds(j,:)>0)
                else
                        ne2=count(redunds(i,:)>0)
                endif
                if (ne2>1) then
                        call gather_web(l, nureads, nnsi)
                else
                        lmark=1
                        allocate(marker(lmark,6))
                        marker(1,:)=recip(k)%v(redunds(l,k),:)
                endif
        endif
        np=np+1
        if (lmark==1) then
                pairs(np,:)=marker(1,:)
        else
                nm6=count(marker(:,6)==maxval(marker(:,6)))
                if (nm6==1) then
                        k=maxloc(marker(:,6),1)
                else
                        k=maxloc(marker(:,5),1,MASK=marker(:,6)==nm6)
                endif
                pairs(np,:)=marker(k,:)
        endif
        do k=1, lmark
                i=marker(k,1); j=marker(k,2)
                redunds(i,:)=0
                redunds(j,:)=0
        enddo
        deallocate(marker)
enddo
print *, '       actual # of pairs = ', np
allocate(AC(nind,np,2), AD(nind,np,2), GC(nind,np), GD(nind,np))
AC=-1
do i=1, np
        AD(:,i,1)=depth(pairs(i,1),:)
        AD(:,i,2)=depth(pairs(i,2),:)
        forall(j=1:nind, AD(j,i,1)>0.and.AD(j,i,2)>0) AC(j,i,1:2)=(/0,1/)
        forall(j=1:nind, AD(j,i,1)>0.and.AD(j,i,2)==0) AC(j,i,1:2)=(/0,0/)
        forall(j=1:nind, AD(j,i,1)==0.and.AD(j,i,2)>0) AC(j,i,1:2)=(/1,1/)
enddo
GC=AC(:,:,1)+AC(:,:,2)
do i=1, np
        forall(j=1:nind, GC(j,i)==-2) GC(j,i)=-1
enddo
GD=AD(:,:,1)+AD(:,:,2)
if (save_opt) then
        open(newunit=u, file='reciprocal_pairs')
        do i=1, np
                write(u,*) pairs(i,:)
        enddo
        close(u)
endif
end subroutine biAssemble2

subroutine gather_web(i, n, m)
implicit none
integer:: i, j, k, l, n, m, n0, n1, inds(2)
integer, allocatable:: v(:)
allocate(v(n))
v=0
v(i)=1
n0=0
n1=sum(v)
do
        if (n0>=n1) exit
        do j=1, n
                if (v(j)==1) then
                        do k=1, m
                                if (redunds(j,k)>0) then
                                        inds=recip(k)%v(redunds(j,k),1:2)
                                        v(inds)=1
                                endif
                        enddo
                        v(j)=2
                endif
        enddo
        n0=n1
        n1=sum(v)
enddo
lmark=0
do j=1, n
        if (v(j)==0) cycle
        lmark=lmark+count(redunds(j,:)>0)
enddo
allocate(marker(lmark,6))
l=0
do j=1, n
        if (v(j)==0) cycle
        do k=1, m
                if (redunds(j,k)==0) cycle
                l=l+1
                marker(l,:)=recip(k)%v(redunds(j,k),:)
        enddo
enddo
end subroutine gather_web

subroutine get_reciprocals(level,list,nl,nc,ksi)
implicit none
integer:: i, j, k, l, level, nt, nr, nl, ksi, nc
integer:: list(nl,nc)
integer, allocatable:: w1(:), w2(:), stat(:)
type(icell), allocatable:: matches(:)
!enumerate reciprocals
nt=count(list(:,3)==level)
if (nt==0) then
        nr=0
        return
endif
allocate(matches(nureads))
matches(:)%n=0
do i=1, nl
        if (list(i,3)/=level) cycle
        j=list(i,1)
        k=list(i,2)
        if (matches(j)%n==0) then
                allocate(matches(j)%v(1))
                matches(j)%v(1)=k
                matches(j)%n=1
        else
                l=matches(j)%n
                allocate(w1(l+1))
                w1(1:l)=matches(j)%v
                w1(l+1)=k
                call move_alloc(w1,matches(j)%v)
                matches(j)%n=l+1
        endif
        if (matches(k)%n==0) then
                allocate(matches(k)%v(1))
                matches(k)%v(1)=j
                matches(k)%n=1
        else
                l=matches(k)%n
                allocate(w1(l+1))
                w1(1:l)=matches(k)%v
                w1(l+1)=j
                call move_alloc(w1,matches(k)%v)
                matches(k)%n=l+1
        endif
enddo
allocate(w1(nureads), stat(nureads))
do i=1, nureads
        if (matches(i)%n==0) cycle
        w1=0
        w1(matches(i)%v)=matches(i)%v
        w1(i)=i
        l=count(w1>0)
        allocate(w2(l))
        w2=pack(w1,w1>0)
        call move_alloc(w2,matches(i)%v)
        matches(i)%n=l
enddo
deallocate(w1)

stat=0
do i=1, nureads
        if (matches(i)%n/=2) cycle
        if (stat(i)/=0) cycle
        j=matches(i)%v(2)
        k=matches(i)%v(1)
        if (k==i.and.matches(j)%n==2) then
                stat(i)=1
                stat(j)=2
        endif
enddo
nr=count(stat==1)
!put them in matrix
allocate(recip(ksi)%v(nr,6))
k=0
do i=1, nureads
        if (stat(i)/=1) cycle
        k=k+1
        j=matches(i)%v(2)
        recip(ksi)%v(k,1:4)=(/i,j,sum(depth(i,1:nind)),sum(depth(j,1:nind))/)
        recip(ksi)%v(k,5)=recip(ksi)%v(k,3)+recip(ksi)%v(k,4)
        recip(ksi)%v(k,6)=min(recip(ksi)%v(k,3),recip(ksi)%v(k,4))
enddo
recip(ksi)%n=nr; recip(ksi)%p=6
deallocate(matches, stat)
end subroutine get_reciprocals

subroutine filtering
! TODO compare GCIII output to Matlab prototype
! TODO remove GCII and GCIV; output of AC/AD
implicit none
integer:: i, j, k, l, u, v, w
integer, allocatable:: minGD(:), w1(:)
real, allocatable:: CR(:), MAF(:), avgGD(:)
logical, allocatable:: MED(:), minGD1(:)
character*16:: fmt1, fmt2
! check for empty columns and monomorphic SNPs in GC
! and take them away
allocate(w1(3))
nm1=0
do i=1, np
        j=count(GC(:,i)>-1)
        if (j==0) cycle ! empty case
        forall(k=1:nind, GC(k,i)>-1) w1(GC(k,i)+1)=w1(GC(k,i)+1)+1
        k=count(w1>0)
        if (k==1) cycle ! monomorphic case
        nm1=nm1+1
        GC(:,nm1)=GC(:,i)
        GD(:,nm1)=GD(:,i)
        AC(:,nm1,1:2)=AC(:,i,1:2)
        AD(:,nm1,1:2)=AD(:,i,1:2)
enddo
deallocate(w1)
print *, '             # of empty markers = ', np-nm1
! compute CR, MAF, minGD, avgGD, MED and minGD1 on selection I
allocate(minGD(nm1), CR(nm1), MAF(nm1), avgGD(nm1), &
        MED(nm1), minGD1(nm1))
MED=.false.
minGD1=.true.
do i=1, nm1
        j=count(GC(:,i)>-1)
        CR(i)=real(j)/real(nind)
        minGD(i)=minval(GD(:,i),MASK=GC(:,i)>-1)
        avgGD(i)=real(sum(GD(:,i),MASK=GC(:,i)>-1))/real(j)
        MAF(i)=.5*real(sum(GC(:,i),MASK=GC(:,i)>-1))/real(j)
        k=get_COVmed(AD(:,i,1),nind)
        l=get_COVmed(AD(:,i,2),nind)
        if (k>=covmed.and.l>=covmed) MED(i)=.true.
        k=count(AD(:,i,1)>0)
        l=minval(AD(:,i,1))
        if (k==1.and.l<covsin) minGD1(i)=.false.
        k=count(AD(:,i,2)>0)
        l=minval(AD(:,i,2))
        if (k==1.and.l<covsin) minGD1(i)=.false.
enddo
! set up next selections (GCII, GCIII, GCIV) and write GCIII to file along with
! CR, MAF, minGD, avgGD
j=count(MED.and.minGD1)
k=count(MED.and.minGD1.and.MAF>=.02.and.CR>=.7)
l=count(MED.and.minGD1.and.MAF>=.02.and.CR==1)
print *, '             # markers (sel II) = ', j
print *, '            # markers (sel III) = ', k
print *, '             # markers (sel IV) = ', l
allocate(GCII(nind,j), GCIII(nind,k,3), GCIV(nind,l))
j=0
k=0
l=0
open(newunit=v, file='sel3.012.markers')
write(v,'(a)') '    sel3#    sel1#       CR      MAF minDepth avgDepth'
do i=1, nm1
        if (MED(i).and.minGD1(i)) then
                j=j+1
                GCII(:,j)=GC(:,i)
                if (MAF(i)<.02.or.CR(i)<.7) cycle
                k=k+1
                GCIII(:,k,1)=GC(:,i)
                GCIII(:,k,2)=AD(:,i,1)
                GCIII(:,k,3)=AD(:,i,2)
                write(v,'(2i9,3f9.3,f9.1)') k, i, CR(i), MAF(i), minGD(i), avgGD(i)
                if (CR(i)<1) cycle
                l=l+1
                GCIV(:,l)=GC(:,i)
        endif
enddo
close(v)
open(newunit=u, file='sel3.012')
open(newunit=v, file='sel3.012.indiv')
open(newunit=w, file='sel3.012.adepth')
write(fmt1,'(a,i12,a)') '(',k,'i3)'
write(fmt2,'(a,i12,a)') '(',k,'i6)'
do i=1, nind
        write(u,fmt1) GCIII(i,:,1)
        write(v,'(a,i3)') indiv(i)%id, indiv(i)%pop
        write(w,fmt2) GCIII(i,:,2)
        write(w,fmt2) GCIII(i,:,3)
enddo
close(u)
close(v)
end subroutine filtering

function get_COVmed(v,n) result(a)
implicit none
integer:: i, j, k, l, n, m, v(n), a
integer, allocatable:: w1(:)
m=count(v>0)
if (m==1) then
        a=v(1)
elseif (m==2) then
        a=minval(v)
else
        l=m/2
        k=maxval(v,MASK=v>0)
        allocate(w1(k))
        w1=0
        forall(i=1:nind, v(i)>0) w1(v(i))=w1(v(i))+1
        i=0
        k=0
        do
                if (i>=l) exit
                k=k+1
                i=i+w1(k)
        enddo
        a=k
        deallocate(w1)
endif
end function get_COVmed

end program
