#!/bin/python

# this file provides auto code generation for mo_util_reorder.F90
import itertools

def code_interface(permutation):
  code =""
  for seq in permutation:
    s = ''.join(seq)
    l = len(seq)
    if (s == ''.join([str(i+1) for i in range(l)])):
      continue
    name = ''.join([str(i+1) for i in range(l)])+'x'+s
    code += '  interface reorder'+name+'\n'
    code += '    module procedure reorder_int_'+name+', reorder_wp_'+name+'\n'
    code += '  end interface'+'\n'
    code += '\n'
  return code

def code_reorder(permutation):
  code =''
  for idx,seq in enumerate(permutation):
    if (idx == 0): continue
    s = ''.join(seq)
    l = len(s)
    idx_left = ['i'+i for i in s]
    idx_right = ['i'+str(i+1) for i in range(l)]
    transpose = int(s[0])>int(s[1])
    idx_right = [i.replace(idx_left[0],':') for i in idx_right]
    idx_right = [i.replace(idx_left[1],':') for i in idx_right]
    idx_left[0] = idx_left[1] = ':'
    func_name = 'reorder_int_'+''.join([str(i+1) for i in range(l)])+'x'+s
    code += '  ! reorder the indecies of 4D array'+'\n'
    code += '  function '+func_name+'(array)'+'\n'
    code += '    integer, dimension('+','.join([':']*l)+'), intent(in) :: array'+'\n'
    code += '    integer, dimension('+','.join(['size(array,dim='+i+')' for i in s])+') :: '+func_name+'\n'
    if (len(idx_left)>2): code += '    integer :: '+','.join(idx_left[2:])+'\n'
    for idx, e in enumerate(reversed(s)):
      if (idx >= (l-2)): continue
      code += '  '*idx+'    do i'+e+' = 1, size(array,dim='+e+')'+'\n'
    if (transpose):
      code += '  '*(l-2)+'    '+func_name+'('+','.join(idx_left)+') = transpose(array('+','.join(idx_right)+'))'+'\n'
    else:
      code += '  '*(l-2)+'    '+func_name+'('+','.join(idx_left)+') = array('+','.join(idx_right)+')'+'\n'
    for idx in reversed(range(l)):
      if (idx >= (l-2)): continue
      code += '  '*idx+'    end do'+'\n'
    code += '  end function'+'\n'
    code += '\n'
    func_name = 'reorder_wp_'+''.join([str(i+1) for i in range(l)])+'x'+s
    code += '  ! reorder the indecies of 4D array'+'\n'
    code += '  function '+func_name+'(array)'+'\n'
    code += '    real(wp), dimension('+','.join([':']*l)+'), intent(in) :: array'+'\n'
    code += '    real(wp), dimension('+','.join(['size(array,dim='+i+')' for i in s])+') :: '+func_name+'\n'
    if (len(idx_left)>2): code += '    integer :: '+','.join(idx_left[2:])+'\n'
    for idx, e in enumerate(reversed(s)):
      if (idx >= (l-2)): continue
      code += '  '*idx+'    do i'+e+' = 1, size(array,dim='+e+')'+'\n'
    if (transpose):
      code += '  '*(l-2)+'    '+func_name+'('+','.join(idx_left)+') = transpose(array('+','.join(idx_right)+'))'+'\n'
    else:
      code += '  '*(l-2)+'    '+func_name+'('+','.join(idx_left)+') = array('+','.join(idx_right)+')'+'\n'
    for idx in reversed(range(l)):
      if (idx >= (l-2)): continue
      code += '  '*idx+'    end do'+'\n'
    code += '  end function'+'\n'
    code += '\n'
  return code

code = """!>
! Module: mo_gas_optics

! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2015,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
! Description: Reorder array indecies
!

module mo_util_reorder

  use mo_rte_kind, only: wp

  implicit none

"""
code += '  ! ----- interface for 2D arrays -----'+'\n'
code += code_interface(itertools.permutations(['1','2']))
code += '  ! ----- interface for 3D arrays -----'+'\n'
code += code_interface(itertools.permutations(['1','2','3']))
code += '  ! ----- interface for 4D arrays -----'+'\n'
code += code_interface(itertools.permutations(['1','2','3','4']))
#code += '  ! ----- interface for 5D arrays -----'+'\n'
#code += code_interface(itertools.permutations(['1','2','3','4','5']))
code += 'contains'+'\n'
code += '\n'
code += '  ! ----- reorder for 2D arrays -----'+'\n'
code += code_reorder(itertools.permutations(['1','2']))
code += '  ! ----- reorder for 3D arrays -----'+'\n'
code += code_reorder(itertools.permutations(['1','2','3']))
code += '  ! ----- reorder for 4D arrays -----'+'\n'
code += code_reorder(itertools.permutations(['1','2','3','4']))
#code += '  ! ----- reorder for 5D arrays -----'+'\n'
#code += code_reorder(itertools.permutations(['1','2','3','4','5']))
code += 'end module'+'\n'
print code
