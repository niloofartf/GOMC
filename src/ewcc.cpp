

Ewald& operator=(Ewald const& rhs){

    for (int i = 0; i < BOXES_WITH_U_NB; i++){

    // These are the vars from Alloc Mem
      kmax[i]               =   rhs.kmax[i];
      imageSize[i]          =   rhs.imageSize[i];
      imageSizeRef[i]       =   rhs.imageSizeRef[i];
      sumRnew[i]            =   rhs.sumRnew[i];
      sumInew[i]            =   rhs.sumInew[i];
      sumRref[i]            =   rhs.sumRref[i];
      sumIref[i]            =   rhs.sumIref[i];
      kx[i]                 =   rhs.kx[i];
      ky[i]                 =   rhs.ky[i];
      kz[i]                 =   rhs.kz[i];
      hsqr[i]               =   rhs.hsqr[i];
      prefact[i]            =   rhs.prefact[i];
      kxRef[i]              =   rhs.kxRef[i];
      kyRef[i]              =   rhs.kyRef[i];
      kzRef[i]              =   rhs.kzRef[i];
      hsqrRef[i]            =   rhs.hsqrRef[i];
      prefactRef[i]         =   rhs.prefactRef[i];

    // The private var
      currentEnergyRecip[i] =   rhs.currentEnergyRecip[i];
    }
      
    imageTotal              =   rhs.imageTotal;
    imageLarge              =   rhs.imageLarge;
      
    alpha                   =   rhs.alpha;
    recip_rcut              =   rhs.recip_rcut;
    recip_rcut_Sq           =   rhs.recip_rcut_Sq;

    particleKind            =   rhs.particleKind;
    particleMol             =   rhs.particleMol;
    particleCharge          =   rhs.particleCharge;

    electrostatic           =   rhs.electrostatic; 
    ewald                   =   rhs.ewald;

    // Since these are external references, I do a shallow copy
    // Note that nothing internally should be changed.
    forcefield              =   rhs.forcefield;
    mols                    =   rhs.mols;
    molLookup               =   rhs.molLookup;

    // I think the references to Box Dim, COM, and syspot can be left alone
    //  They should still point to the old address, which is copy constructed from the replica
    currentCoords           =   rhs.currentCoords;
    currentCOM              =   rhs.currentCOM;
    systemPotRef            =   rhs.sysPotRef;

}


