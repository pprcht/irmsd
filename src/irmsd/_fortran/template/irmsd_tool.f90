subroutine irmsd_tool(fname1,fname2,iinversion)                        
!*******************************************************               
!* irmsd_tool                                                          
!* Standalone implementation to compare two structures                 
!* with the iRMSD method.                                              
!* This implementation should be called only on its own,               
!* for ensemble-based processing see the CREGEN file                   
!*******************************************************               
  use crest_parameters                                                 
  use strucrd                                                          
  use axis_module                                                      
  use irmsd_module                                                     
  use canonical_mod                                                    
  implicit none                                                        
  character(len=*),intent(in) :: fname1                                
  character(len=*),intent(in) :: fname2                                
  integer,intent(in) :: iinversion                                     
  type(coord) :: mol,ref                                               
  real(wp) :: rmsdval,tmpd(3),tmpdist                                  
  integer :: i,ich                                                     
  type(rmsd_cache) :: rcache                                           
  type(canonical_sorter) :: canmol                                     
  type(canonical_sorter) :: canref                                     
  logical :: mirror                                                    
  logical,parameter :: debug = .false.                                 
                                                                       
  write (stdout,*) 'iRMSD algorithm'                                   
  write (stdout,*) 'reference: ',trim(fname1)                          
  write (stdout,*) 'processed: ',trim(fname2)                          
  write (stdout,*)                                                     
                                                                       
  !> read the geometries                                               
  call ref%open(trim(fname1))                                          
  call mol%open(trim(fname2))                                          
                                                                       
  !> move ref to CMA and align rotational axes                         
  call axis(ref%nat,ref%at,ref%xyz)                                    
                                                                       
  !> allocate memory                                                   
  call rcache%allocate(ref%nat)                                        
                                                                       
  !> canonical atom ranks                                              
  call canref%init(ref,invtype='apsp+',heavy=.false.)                  
  !call canref%add_h_ranks(ref)                                        
  rcache%stereocheck = .not. (canref%hasstereo(ref))                   
  call canref%shrink()                                                 
  write (stdout,*) 'false enantiomers possible?: ',rcache%stereocheck  
  select case (iinversion)                                             
  case (0)                                                             
    mirror = .true.                                                    
  case (1)                                                             
    mirror = .true.                                                    
    rcache%stereocheck = .true.                                        
  case (2)                                                             
    mirror = .false.                                                   
    rcache%stereocheck = .false.                                       
  end select                                                           
  write (stdout,*) 'allow inversion?:            ',mirror              
                                                                       
  call canmol%init(mol,invtype='apsp+',heavy=.false.)                  
  !call canmol%add_h_ranks(mol)                                        
  call canmol%shrink()                                                 
                                                                       
  !> check if we can work with the determined ranks                    
  if (checkranks(ref%nat,canref%rank,canmol%rank)) then                
    write (stdout,*) 'using canonical atom identities as rank backend' 
    rcache%rank(:,1) = canref%rank(:)                                  
    rcache%rank(:,2) = canmol%rank(:)                                  
    if (debug) then                                                    
      write (*,*) 'iRMSD ranks:'                                       
      write (*,*) 'atom',' rank('//fname1//')',' rank('//fname2//')'   
      do i = 1,ref%nat                                                 
        write (*,*) i,rcache%rank(i,1),rcache%rank(i,2)                
      end do                                                           
      write (*,*)                                                      
    end if                                                             
  else                                                                 
    !> if not, fall back to atom types                                 
    write (stdout,*) 'using atom types as rank backend'                
    call fallbackranks(ref,mol,ref%nat,rcache%rank)                    
  end if                                                               
                                                                       
  call min_rmsd(ref,mol,rcache=rcache,rmsdout=rmsdval,align=.true.)    
                                                                       
  !> write the rotated and shifted coordinates to one file             
  open (newunit=ich,file='irmsd.xyz')                                  
  call ref%append(ich)                                                 
  call mol%append(ich)                                                 
  close (ich)                                                          
  write (stdout,*)                                                     
  write (stdout,*) 'aligned structures written to irmsd.xyz'           
  write (stdout,*)                                                     
                                                                       
  do i = 1,mol%nat                                                     
    tmpd(:) = (mol%xyz(:,i)-ref%xyz(:,i))**2                           
    tmpdist = sqrt(sum(tmpd(:)))*autoaa                                
    if (tmpdist > 0.01_wp) then                                        
      write (*,*) i,mol%at(i),tmpdist                                  
    end if                                                             
  end do                                                               
                                                                       
  rmsdval = rmsdval*autoaa                                             
  write (*,'(1x,a,f16.8)') 'Calculated iRMSD (Å):',rmsdval             
                                                                       
  return                                                               
end subroutine irmsd_tool                                              

