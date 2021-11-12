remark topology file for dipolar coupling calculations
set message off echo off end

autogenerate 
  angles=true
  dihedrals=false
end

mass XX   100
mass YY	  100
mass ZZ   100
mass OO   100

residue ANI
  group
   atom X  type=XX  charge=0.0  end
   atom Y  type=YY  charge=0.0  end
   atom Z  type=ZZ  charge=0.0  end
   atom OO type=OO  charge=0.0  end

   bond X OO
   bond Y OO
   bond Z OO
end

residue XAN
  group
   atom X  type=XX  charge=0.0  end
   atom Y  type=YY  charge=0.0  end
   atom Z  type=ZZ  charge=0.0  end
   atom OO type=OO  charge=0.0  end

   bond X OO
   bond Y OO
   bond Z OO
end

residue DAN
  group
   atom X  type=XX  charge=0.0  end
   atom Y  type=YY  charge=0.0  end
   atom Z  type=ZZ  charge=0.0  end
   atom OO type=OO  charge=0.0  end

   bond X OO
   bond Y OO
   bond Z OO
end

set message on echo on end
