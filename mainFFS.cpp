/*********************************************************************

A liquid argon drop spreading on structurelss patterned substrate 
  

**********************************************************************/


#include "headfileFFS.h"
double const pi = 4*std::atan(1.0);

 

int main ()
{
   
 //int starttime, stoptime, timeused;
 //starttime = clock();


 nMol = 4*VProd (initUcell);//fcc lattice
 HalfNmol = 12433;
 
 AllocArrays ();
 DropInitReading();// to read initial drop information
 //HemisphereInit();//to cut off the sphere
 //PrinthalfdropletCoordInit();//initial hemisphere

 //contact_angle[0]=ContactAngleCalculation(); //initial contact angle;
 //PrintContactAngle(contact_angle[0]);
 //cout<<"initial contact angle="<<contact_angle[0]<<endl;
 //PrintInterface(); //output of discrete interface
 //exit(0);
 


 
 dropspreading();
 
 PrintEquilCoordVelAcc();//information at equilibrium state
 PrintVelocityDistribution();
 //stoptime = clock();
 //timeused =stoptime - starttime;
 //timecounting(timeused);

 return (0);
}



 
void AllocArrays ()
{

  AllocMem (mol, nMol, Mol);
  AllocMem (cellList, VProd (cells) + nMol, int);
  AllocMem (histVel,sizeHistVel,double);//velocity distribution
  AllocMem (halfmol, HalfNmol, Mol);
 // AllocMem (contact_angle,int(stepLimit/stepOutPut),double);
  AllocMem(Fitatoms, nMol, Mol);
  AllocMem(CentralDroplet, nMol, Mol);
  AllocMem(LayerAtoms, nMol, Mol);
  AllocMem(fitmol, nMol, Mol);
  AllocMem(Uppmol, nMol, Mol);
  AllocMem (halfmolGap, nMol, Mol);
			

 
}

/*********---to read initial coordinates, velocites, or acceleration------------------****/
void DropInitReading()
{

    int n=0;
	ifstream inCoord("Dropinitcoord.txt",ios::in|ios::binary);
	ifstream inVel("Dropinitvel.txt",ios::in|ios::binary);
	ifstream inAcc("Dropinitacc.txt",ios::in|ios::binary);



	if((!inCoord)||(!inVel)||(!inAcc))
	{
		cout<<"could not open file"<<endl;
	}

	for (n=0;n<HalfNmol;n++)
	{
		inCoord>>halfmol[n].r.x;
		inCoord>>halfmol[n].r.y;
		inCoord>>halfmol[n].r.z;

		inVel>>halfmol[n].rv.x;
		inVel>>halfmol[n].rv.y;
		inVel>>halfmol[n].rv.z;

		inAcc>>halfmol[n].ra.x;
		inAcc>>halfmol[n].ra.y;
		inAcc>>halfmol[n].ra.z;
	}
	inCoord.close();
	inVel.close();
	inAcc.close();


}




/*********-----------------to get hemisphere by cutting----------------------------**********/

void HemisphereInit()
{

	int n;
	int m=0;
	DO_MOL {
		if(mol[n].r.z>=0.0&& mol[n].r.z<=2.*Radius && fabs(mol[n].r.y)<=2.*Radius)
		{
		
			halfmol[m].r.x=mol[n].r.x;
			halfmol[m].r.y=mol[n].r.y;
			halfmol[m].r.z=mol[n].r.z;

			halfmol[m].rv.x=mol[n].rv.x;
			halfmol[m].rv.y=mol[n].rv.y;
			halfmol[m].rv.z=mol[n].rv.z;

			halfmol[m].ra.x=mol[n].ra.x;
			halfmol[m].ra.y=mol[n].ra.y;
			halfmol[m].ra.z=mol[n].ra.z;

			m=m+1;
		}
		
	};


  HalfNmol=m;

}



void HemisphereInitContinue()
{
	int n;
	int m=0;
	//halfmol=mol;
	ifstream inCoord("Dropinitcoord1.txt",ios::in|ios::binary);

	for (n=0;n<HalfNmol;n++)
	{
		inCoord>>halfmol[n].r.x;
		inCoord>>halfmol[n].r.y;
		inCoord>>halfmol[n].r.z;

	
	}
	inCoord.close();

}





/*void PrinthalfdropletCoordInit()
{

	int n;
	
	ofstream outfile_halfdroplet("halfdropletinit.txt",ios::out);
	

	outfile_halfdroplet<<"HalfNmol="<<HalfNmol<<endl;
	outfile_halfdroplet<<"r.x"<<setw(14)<<"r.y"<<setw(14)<<"r.z"<<endl;
	for (n=0;n<HalfNmol;n++)
	{
		outfile_halfdroplet<<halfmol[n].r.x<<setw(14)<<halfmol[n].r.y<<setw(14)<<halfmol[n].r.z<<endl;
	 }


	outfile_halfdroplet<<"rv.x"<<setw(14)<<"rv.y"<<setw(14)<<"rv.z"<<endl;
	for (n=0;n<HalfNmol;n++)
	{
		outfile_halfdroplet<<halfmol[n].rv.x<<setw(14)<<halfmol[n].rv.y<<setw(14)<<halfmol[n].rv.z<<endl;
	  
	 }

	outfile_halfdroplet<<"ra.x"<<setw(14)<<"ra.y"<<setw(14)<<"ra.z"<<endl;
	for (n=0;n<HalfNmol;n++)
	{
		outfile_halfdroplet<<halfmol[n].ra.x<<setw(14)<<halfmol[n].ra.y<<setw(14)<<halfmol[n].ra.z<<endl;
	  
	 }
    
	outfile_halfdroplet.close();
}
*/



void PrinthalfdropletCoordMove()
{

	int n;
	
	ofstream outfile_halfdroplet("halfdropletMove.txt",ios::out);
	

	outfile_halfdroplet<<"HalfNmol="<<HalfNmol<<endl;
	outfile_halfdroplet<<"r.x"<<setw(14)<<"r.y"<<setw(14)<<"r.z"<<endl;
	for (n=0;n<HalfNmol;n++)
	{
		outfile_halfdroplet<<halfmol[n].r.x<<setw(14)<<halfmol[n].r.y<<setw(14)<<halfmol[n].r.z<<endl;
	 }
    
	outfile_halfdroplet.close();
}








 

//*************************structureless simulation******************/



void dropspreading()
{

  //Structureless_SystemMove();//to move hemispherical drop up InitialMove by considering no-penetration
  density=1.0; //initial density
  VSCopy (region, 1. / pow (density/4., 1./3.), initUcell);//rou*L1*L2*L3=nMol
  VSCopy (simulationregion, powl(2.0,1./3), region);//to let lattice fill half the volume of a cube

   
  velMag = sqrt(NDIM * (1. - 1. / HalfNmol) * temperature);//velMag*velMag=d*T
  
  stepCount = 0;
   srand ( time(NULL) );
	int iSecret;
	iSecret = rand() % 100 + 1;
	randSeedP=iSecret;
  HalfdropInitVels();
  HalfdropInitAccels();

  PrinthalfdropletCoordInit();

   
 
  simulationregion.y=NSTEPS*StepWidth+NSTEPS*StepGap;	
  
  if(StepHeight==0) {simulationregion.z=0.5*simulationregion.y+InitialMove;} 
  else { simulationregion.z=simulationregion.y+2*StepHeight+InitialMove;}

  VSCopy (halfdrop_cells, 1. / rCut, simulationregion);  //(v2).x = (s1) * (v1).x
  AllocMem (halfdrop_cellList, VProd (halfdrop_cells) + HalfNmol, int);
  
  //stepLimit=800000; 
   
  moreCycles = 1;
  while (moreCycles) 
  { 
	// cout<<"stepCount="<<stepCount<<endl;
	Structureless_SingleStep (); 	
	 CountDensity();
	 CenterOfMassCalc();
	

     if (stepCount >= stepLimit) moreCycles = 0;
	// if (stepCount%(250*stepOutPut)== 0)
	 
	
			

		 //HalfdropInitVels();
		 // PrintStructurelessCoordEvalution();

	//if (stepCount== 46671)
	// {
	//	PrintStructurelessCoordEvalution();
	//	PrintDensity();
	//}

	/* if (stepCount%(200*stepOutPut)== 0)
	{
		
		
		PrintTotalEnergy();
		
	 }*/

	 //if (stepCount%stepOutPut== 0)
	// {
		//NCA=NCA+1;
		//PrintStructurelessCoordEvalution();
		//PrintTotalEnergy();
		//forceztotal();
		//cout<<stepCount<<setw(14)<<totalU<<endl;
		//contact_angle[NCA]=ContactAngleCalculation();
		//PrintContactAngle(contact_angle[NCA]);
		//cout<<"stepCount="<<stepCount<<endl;
		//cout<<"contact_angle[NCA]="<<contact_angle[NCA]<<endl;
		//cout<<"CA difference="<<fabs(contact_angle[NCA]-contact_angle[NCA-1])<<endl;
		//if(fabs(contact_angle[NCA]-contact_angle[NCA-1])<0.01){moreCycles=0;} //using contact angle criteria to terminate the running
		//if(fabs(contact_angle[NCA])<3.0) {moreCycles=0;}
	// }
	 
	  
  }


}




void CountDensity()
{
	int nGap1=0;
	int n; 
	int m=0;
	for (n=0;n<HalfNmol;n++)
	{
		if (((halfmol[n].r.y)<= (StepGap/2)) && ((halfmol[n].r.y) >= (-StepGap/2)) && ((halfmol[n].r.z) <= (StepHeight))   )


		{
			nGap1=nGap1+1;
			halfmolGap[m].r.x=halfmol[n].r.x;
			halfmolGap[m].r.y=halfmol[n].r.y;
			halfmolGap[m].r.z=halfmol[n].r.z;

			halfmolGap[m].rv.x=halfmol[n].rv.x;
			halfmolGap[m].rv.y=halfmol[n].rv.y;
			halfmolGap[m].rv.z=halfmol[n].rv.z;

			halfmolGap[m].ra.x=halfmol[n].ra.x;
			halfmolGap[m].ra.y=halfmol[n].ra.y;
			halfmolGap[m].ra.z=halfmol[n].ra.z;

			m=m+1;
					
		}

	}
	DensityGap1=nGap1/(4*StepGap*StepHeight);
	
	double initDensity=0.34;
	double DesiredDensity=1.65;
	double DesiredDensity2=1.6;
	double DesiredDensity3=1.45;
	double DesiredDensity4=1.4;
	double DesiredDensity5=1.35;
	double DesiredDensity6=1.3;

	



	
	RoundDensity=roundMaker(DensityGap1);
	PrintDensity();
	
	if ((RoundDensity < initDensity))
	{		PrintFinish();
	        exit(0);

	}



  if ((RoundDensity >= DesiredDensity4) & (stepDesired4 > stepCount) )

	{
		stepDesired4=stepCount;
		PrintStructurelessCoordEvalution4();
		
	}
/*
    if ((RoundDensity >= DesiredDensity5) & (stepDesired5 > stepCount) )

	{
		stepDesired5=stepCount;
		PrintStructurelessCoordEvalution5();
		
	}

	  if ((RoundDensity >= DesiredDensity6) & (stepDesired6 > stepCount) )

	{
		stepDesired6=stepCount;
		PrintStructurelessCoordEvalution6();
		
	}

*/


  if ((RoundDensity >= DesiredDensity3) & (stepDesired3 > stepCount) )

	{
		stepDesired3=stepCount;
		PrintStructurelessCoordEvalution3();
		cout<<stepDesired3<<endl;
		

	}


	if ((RoundDensity >= DesiredDensity2) & (stepDesired2 > stepCount) )

	{
		stepDesired2=stepCount;
		PrintStructurelessCoordEvalution2();
		cout<<stepDesired2<<endl;
		

	}


	if (RoundDensity >= DesiredDensity)
	{
		PrintStructurelessCoordEvalution();
		exit(0);

	}


	
	
}

void PrintFinish()
{
	
	ofstream outfile_Finish("Finish.txt",ios::app);
	outfile_Finish<<stepCount*deltaT<<setw(14)<<DensityGap1<<setw(14)<<RoundDensity<<endl;
	outfile_Finish.close();

}


double roundMaker(double DensityGap1)
{
	double DeltaValue;
	DeltaValue=DensityGap1*100-floor(DensityGap1*100);
	//cout<<DeltaValue<<endl;
	double RoundValue;
	if (DeltaValue >= 0.5 )
	{
		RoundValue=DensityGap1-(DeltaValue/100)+0.01;
		//cout<<RoundValue<<endl;

	}
	else
	{
		RoundValue=DensityGap1-DeltaValue/100;
	}


return(RoundValue);
}

void PrintDensity()
{
	
	ofstream outfile_DensityGap("DensityGap.txt",ios::app);
	outfile_DensityGap<<stepCount*deltaT<<setw(14)<<DensityGap1<<setw(14)<<RoundDensity<<endl;
	outfile_DensityGap.close();

}
 



void PrintGapCoord(int nGap1)
{
    int n;
	ofstream outfile_dropevolutionG("spreading_processG.txt",ios::app);
	ofstream outfile_velocityG("velocityG.txt",ios::app);
	ofstream outfile_accelerationG("accelerationG.txt",ios::app);

	outfile_dropevolutionG<<"time="<<setw(14)<<stepCount*deltaT<<endl;
	outfile_dropevolutionG<<"x"<<setw(14)<<"y"<<setw(14)<<"z"<<endl;

    outfile_velocityG<<"time="<<setw(14)<<stepCount*deltaT<<endl;
	outfile_velocityG<<"vx"<<setw(14)<<"vy"<<setw(14)<<"vz"<<endl;

	 outfile_accelerationG<<"time="<<setw(14)<<stepCount*deltaT<<endl;
	outfile_accelerationG<<"ax"<<setw(14)<<"ay"<<setw(14)<<"az"<<endl;




	//outfile_halfdroplet<<"x"<<setw(14)<<"y"<<setw(14)<<"z"<<setw(14)<<"rv.x"<<setw(14)<<"rv.y"<<setw(14)<<"rv.z"<<setw(14)<<"ra.x"<<setw(14)<<"ra.y"<<setw(14)<<"ra.z"<<endl;	 	  
	for (n=0;n<nGap1;n++)
	{
		outfile_dropevolutionG<<halfmolGap[n].r.x<<setw(14)<<halfmolGap[n].r.y<<setw(14)<<halfmolGap[n].r.z<<endl;
		outfile_velocityG<<halfmolGap[n].rv.x<<setw(14)<<halfmolGap[n].rv.y<<setw(14)<<halfmolGap[n].rv.z<<endl;
		outfile_accelerationG<<halfmolGap[n].ra.x<<setw(14)<<halfmolGap[n].ra.y<<setw(14)<<halfmolGap[n].ra.z<<endl;
	 }
	outfile_dropevolutionG.close();
	outfile_velocityG.close();
	outfile_accelerationG.close();


}









/*
void Structureless_SingleStep ()//checking here
{

  ++ stepCount;
  timeNow = stepCount * deltaT;

 Structureless_LeapfrogStep (1);
 Structureless_ApplyBoundaryCond();
 Structureless_ComputeForcesCellDiv();
 Structureless_LeapfrogStep(2); 
 if (stepCount%(10*stepOutPut)== 0) VelocityRescaling();
 
 }
 */
 

void Structureless_SingleStep ()//checking here
{

  ++ stepCount;
  timeNow = stepCount * deltaT;
  PredictorStep ();
  Structureless_ApplyBoundaryCond();
  Structureless_ComputeForcesCellDiv();
   NoseHooverThermostat();
  CorrectorStep ();
  
  Structureless_ApplyBoundaryCond();
 
  
}

/*
void Structureless_SingleStep ()//checking here
{

  ++ stepCount;
  timeNow = stepCount * deltaT;

 velocity_Verlet (1);
 Structureless_ApplyBoundaryCond();
 Structureless_ComputeForcesCellDiv();
 velocity_Verlet(2); 
 
 
 }

 */


void Structureless_LeapfrogStep (int part)
{
  int n;

  if (part == 1) 
  {
	 for (n=0;n<HalfNmol;n++)
	 {
		VVSAdd (halfmol[n].rv, 0.5*deltaT, halfmol[n].ra);
		VVSAdd (halfmol[n].r, deltaT, halfmol[n].rv);
     }
  } 

  else
  {
     for (n=0;n<HalfNmol;n++)
	 {
		VVSAdd (halfmol[n].rv, 0.5*deltaT, halfmol[n].ra);
	 }
  }
}

///////

void velocity_Verlet(int part)
{
  int n;
  vvSum=0;
  double tempa;
  double nu=0.1;
  int n_dim=3;
   double   mue=0;
  if (part == 1) 
  {
	 for (n=0;n<HalfNmol;n++)
	 {
		VVSAdd (halfmol[n].r, deltaT, halfmol[n].rv);
		VVSAdd (halfmol[n].r, 0.5*deltaT*deltaT, halfmol[n].ra);
		VVSAdd (halfmol[n].rv, 0.5*deltaT, halfmol[n].ra);
		
     }
  } 

  else
  {
     for (n=0;n<HalfNmol;n++)
	 {
		VVSAdd (halfmol[n].rv, 0.5*deltaT, halfmol[n].ra);
		vvSum += VLenSq (halfmol[n].rv);
		
	 }
	 tempa=vvSum/(n_dim*(HalfNmol-1));
	 double  sigma=sqrt(tempa);
	srand ( time(NULL) );
	int iSecret;
	iSecret = rand() % 100 + 1;
	 //double probab=gasdev(sigma,mue);
	double probab=RandR();
	//cout<<deltaT*nu<<setw(14)<<probab<<endl;
	 if ((deltaT*nu)>probab)
	 {

		 HalfdropInitVels ();
		 
		 
		 

	 }


  }
	   
}
/////






void Structureless_ApplyBoundaryCond()
{

	int n;
	for (n=0;n<HalfNmol;n++)
	{
		
		VWrap(halfmol[n].r,x);
		VWrap(halfmol[n].r,y);
 	 }
	 
}

void Structureless_ComputeForcesCellDiv()
{
	VecR dr,invWid,rs,shift;
	VecI cc,m1v,m2v,vOff[]=OFFSET_VALS;
	double fcVal,rr,rrCut,rri,rri3;
	double dur=0.0,dua=0.0;
	int c,j1,j2,m1,m1x,m1y,m1z,m2,n,offset;
	


	int nstep;
	double yl, ys,zl,zs;
	double term1,term2,term3,term4,term5,term6,term7,term8;


	rrCut=Sqr(rCut);

    for (n=0;n<HalfNmol;n++)
	{
	 VZero(halfmol[n].ra);
	}
	

	uSum=0.0;//initialization
	virSum=0.0;//initialization
	VDiv(invWid,halfdrop_cells,simulationregion);//physical dimensions of each cells (v1).x = (v2).x / (v3).x, etc. 
	for(n=HalfNmol;n<HalfNmol+VProd(halfdrop_cells);n++){ halfdrop_cellList[n]=-1;}//initialized a pointer associated with cells

  	for (n=0;n<HalfNmol;n++)
	{
	
		 //VSAdd(rs,halfmol[n].r,0.5,simulationregion); //(v1).x = (v2).x + (s3) * (v3).x, etc.
		rs.x=halfmol[n].r.x+0.5*simulationregion.x; //to shift coordinate system to the bottom of the simulation region
		rs.y=halfmol[n].r.y+0.5*simulationregion.y;
		rs.z=halfmol[n].r.z;

	
		VMul(cc,rs,invWid);//cc=rs*invWid::rescaled by invWid and deterine the cell position for each atom
		c=VLinear(cc,halfdrop_cells)+HalfNmol;// VLinear(p, s) =((p).z * (s).y + (p).y) * (s).x + (p).x) //
	  	halfdrop_cellList[n]=halfdrop_cellList[c];//different atoms belonging to the same cell
		halfdrop_cellList[c]=n;//one per cell, point to the first atom in each cell
	}

   

	for(m1z=0;m1z<halfdrop_cells.z;m1z++)//loop on cells-z component
	{

		for(m1y=0;m1y<halfdrop_cells.y;m1y++)//loop on cells-y component
		{
			for(m1x=0;m1x<halfdrop_cells.x;m1x++)//loop on cells-x component
			{
				VSet(m1v,m1x,m1y,m1z);//using a vector to characterize a cell  (m1v).x =m1x,etc.  
				m1=VLinear(m1v,halfdrop_cells)+HalfNmol;//the atoms  and the numbering cells
				for(offset=0;offset<N_OFFSET;offset++)//(check 14 neighbor cells)
				{
					VAdd(m2v,m1v,vOff[offset]);//adding current cells to checking neighbor cell
					VZero(shift);//initial shift={0,0,0}
					HalfVCellWrap(x);
					HalfVCellWrap(y);
					if(m2v.z<0||m2v.z>=halfdrop_cells.z) continue;
                 	m2=VLinear(m2v,halfdrop_cells)+HalfNmol;//to find the number of atoms in all neighbor cells
					 
                	for (j1 = halfdrop_cellList[m1]; j1 >= 0; j1= halfdrop_cellList[j1])
					{
					

						for (j2 = halfdrop_cellList[m2]; j2 >= 0; j2= halfdrop_cellList[j2])
						{
							
							if(m1!=m2||j2<j1)
							{
									
								VSub(dr,halfmol[j1].r,halfmol[j2].r);//dr.x=mol[j1].r.x-mol[j2].r.x;
								VVSub(dr,shift);//#define VVSub(v1, v2)  VSub (v1, v1, v2)
							 	rr=VLenSq(dr);
								if(rr<rrCut)
								{
									rri=1/rr;
									rri3=Cube(rri);
									fcVal=48.*rri3*(rri3-0.5)*rri;
									VVSAdd(halfmol[j1].ra,fcVal,dr);
									VVSAdd(halfmol[j2].ra,-fcVal,dr);
									
								}
							}		
						}
					}
				}
			}
		}
	}

    for (n=0;n<HalfNmol;n++)
	{

	
		dr.x=0.0;
		dr.y=0.0;
		//dr.z=halfmol[n].r.z-0.5*simulationregion.z;
		dr.z=halfmol[n].r.z-simulationregion.z; //changed
		fcVal=-4./5*M_PI*Egs*rous*Lgs12/powl(dr.z,10.0);//repulsive upper wall 
		//if(fcVal>10000) {cout<<"here1 from top wall";}
		halfmol[n].ra.z=halfmol[n].ra.z+fcVal; 

	
	}




   //to calcualte the interaction between a liquid atom and a grooved surface
	for (n=0;n<HalfNmol;n++)
	{

		dr.x=0.0;
		dr.y=0.0;
	
       /*******************************----------------------****/
		//**this part is to consider the step walls' effects***//

		for(nstep=-hstepGhost;nstep<NSTEPS+hstepGhost;nstep++) //practical steps+ghost steps
		{
		
			if(StepHeight<=0) {break;} //in case of flat surface
			ys=-NSTEPS/2.*(StepWidth+StepGap)+0.5*StepGap-halfmol[n].r.y+nstep*(StepGap+StepWidth);
			yl=ys+StepWidth;
			zs=-halfmol[n].r.z;
			zl=zs+StepHeight;
			if(fabs(yl)>=tol && fabs(ys)>=tol&& fabs(zs)>=tol && fabs(zl)>=tol)
			{
				
				
				//to calcualte force of y component
				

				if(fabs(yl)>DCut) //atom away from the step side of yl
				{
					term1=0.0;
					term3=0.0;
					term5=0.0;
					term7=0.0;
				}
				else
				{
					term1=FUrep(yl,zl);
					term3=FUrep(yl,zs);
					term5=FUatt(yl,zl);
					term7=FUatt(yl,zs);
				}
				if(fabs(ys)>DCut) //steele potential  //atom away from the step side of ys
				{
					term2=0.0;
					term4=0.0;
					term6=0.0;
					term8=0.0;
				
				}
				else
				{

					term2=FUrep(ys,zs);
					term4=FUrep(ys,zl);
					term6=FUatt(ys,zs);
					term8=FUatt(ys,zl);
				
				}

				/*term1=FUrep(yl,zl);
				term3=FUrep(yl,zs);
				term5=FUatt(yl,zl);
				term7=FUatt(yl,zs);
				term2=FUrep(ys,zs);
				term4=FUrep(ys,zl);
				term6=FUatt(ys,zs);
				term8=FUatt(ys,zl);*/

				dur=term1+term2-term3-term4;
				dua=term5+term6-term7-term8;
				fcVal=4.0*Egs*rous*(Lgs12*dur-Lgs6*dua); //sign here
				//if(fcVal>10000) {cout<<"y nonsingular";}
				halfmol[n].ra.y=halfmol[n].ra.y+fcVal;  // to update y component force 



               //to calculate z-component force 
				if(fabs(zl)>DCut) //atom away from the step side of yl
				{
					term1=0.0;
					term3=0.0;
					term5=0.0;
					term7=0.0;
				}
				else
				{
					term1=FUrep(zl,yl);
					term3=FUrep(zl,ys);
					term5=FUatt(zl,yl);
					term7=FUatt(zl,ys);
				}
				if(fabs(zs)>DCut) //steele potential  //atom away from the step side of ys
				{
					term2=0.0;
					term4=0.0;
					term6=0.0;
					term8=0.0;
				
				}
				else
				{

					term2=FUrep(zs,ys);
					term4=FUrep(zs,yl);
					term6=FUatt(zs,ys);
					term8=FUatt(zs,yl);
				
				}
				/*term1=FUrep(zl,yl);
				term3=FUrep(zl,ys);
				term5=FUatt(zl,yl);
				term7=FUatt(zl,ys);
				term2=FUrep(zs,ys);
				term4=FUrep(zs,yl);
				term6=FUatt(zs,ys);
				term8=FUatt(zs,yl);*/

				dur=term1+term2-term3-term4;
				dua=term5+term6-term7-term8;
				fcVal=4.0*Egs*rous*(Lgs12*dur-Lgs6*dua);
				//if(fcVal>10000) {cout<<"z nonsingular";}
				halfmol[n].ra.z=halfmol[n].ra.z+fcVal;  // to update z component force 
			}
			else
			{
				halfmol[n].ra.y=halfmol[n].ra.y+AForce(yl,ys,zl,zs);  // to update y component force
			//	if(fcVal>10000) {cout<<"y singular";}
				halfmol[n].ra.z=halfmol[n].ra.z+AForce(zl,zs,yl,ys); // to update z component force
			//	if(fcVal>10000) {cout<<"z singular";}
			}

		}  //end of step loop




    	//Flat base below the steps	
		dr.z=halfmol[n].r.z;
		//if(fabs(dr.z)>=DCut){continue;}
		dua=1./(2*powl(dr.z,4.0));
		dur=1./(5*powl(dr.z,10.0));

		fcVal=4.0*M_PI*Egs*rous*(Lgs12*dur-Lgs6*dua);
	//	cout<<"fcVal="<<fcVal<<"z with flat"<<endl;
	//	if(fcVal>10000) {cout<<"z with flat";}
	    halfmol[n].ra.z=halfmol[n].ra.z+fcVal;
	}
}



  
  

double Sign(double y)
{
	if(y<0.0) return (-1.0);
	else return (1.0);
}

void ApplyThermostat ()
{
  VecR vt;
  double  s1, s2, vFac;
  int n;

  s1 =0;
  s2 = 0.;
  for (n=0;n<HalfNmol;n++) {
    VSAdd (vt, halfmol[n].rv, 0.5 * deltaT, halfmol[n].ra);
    s1 += VDot (vt, halfmol[n].ra);
    s2 += VLenSq (vt);
  }
  vFac = - s1 / s2;
  for (n=0;n<HalfNmol;n++){
    VSAdd (vt, halfmol[n].rv, 0.5 * deltaT,halfmol[n].ra);
    VVSAdd (halfmol[n].ra, vFac, vt);
  }
}



void ApplyThermostat2()
{
  
  double  s1, s2, vFac;
  int n;

  s1 =0;
  s2 = 0.;
  for (n=0;n<HalfNmol;n++) {
    s1 += VDot (halfmol[n].rv, halfmol[n].ra);
    s2 += VLenSq (halfmol[n].rv);
  }
  vFac = - s1 / s2;
  for (n=0;n<HalfNmol;n++){
    VVSAdd (halfmol[n].ra, vFac, halfmol[n].rv);
  }
}

void CenterOfMassCalc()
{
	//VecR SumCenter;
	VZero (SumCenter);
	int n;
	for (n=0;n<HalfNmol;n++)
	{
		VVAdd (SumCenter, halfmol[n].r);

	}
    ZCenterofMass=SumCenter.z/HalfNmol;
	RoundCenter=roundMaker(ZCenterofMass);
	//cout<<RoundCenter<<endl;

	PrintCenterofMass();

}


void PrintCenterofMass()
{
	ofstream outfile_COM("COM.txt",ios::app);
	outfile_COM<<stepCount*deltaT<<setw(14)<<ZCenterofMass<<setw(14)<<RoundCenter<<setw(14)<<endl;
	outfile_COM.close();
}


void NoseHooverThermostat()
{

  double s2;
  int n;

  ///////////////////substract the center of mass
   VZero (vSum);
	for (n=0;n<HalfNmol;n++)
	{
		VVAdd (vSum, halfmol[n].rv);
	}
	for (n=0;n<HalfNmol;n++) { VVSAdd (halfmol[n].rv, - 1. / HalfNmol, vSum);}
 /////////////////////////////////////////////////////////////////////////////////

  s2 = 0.;
  for (n=0;n<HalfNmol;n++) {
       s2 += VLenSq (halfmol[n].rv);
  }
  gama=gama+deltaT*(-3*HalfNmol*temperature+s2)/Q;
  for (n=0;n<HalfNmol;n++){
    VVSAdd (halfmol[n].ra, -gama, halfmol[n].rv);
  }
  

}

void VelocityRescaling()
{
	int n;
	double lamda;
	
	double uk=0;
	for (n=0;n<HalfNmol;n++)
	{
		 uk=uk+VLenSq(halfmol[n].rv);		
     }
	lamda=sqrtl(temperature*3*HalfNmol/uk);
	for (n=0;n<HalfNmol;n++)
	{
		VScale (halfmol[n].rv, lamda);	
     }
	//cout<<"lamda="<<lamda<<endl;
	//if(stepCount>=10000) exit(0);

}


void AdjustTemp()
{
	double vFac;
	int n;
	vvSum=0;
	for (n=0;n<HalfNmol;n++) 
	{
       vvSum += VLenSq (halfmol[n].rv);


	}
	vFac=velMag/sqrt(vvSum/HalfNmol);
	for (n=0;n<HalfNmol;n++) 
	{
       VScale(halfmol[n].rv,vFac);
	}


}

void AdjustVelocity()
{
	int n=0;
	VZero (vSum);
	for (n=0;n<HalfNmol;n++) 
	{
	
		VVAdd (vSum, halfmol[n].rv);
	}
	for (n=0;n<HalfNmol;n++) { VVSAdd (halfmol[n].rv, - 1. / HalfNmol, vSum);} //to keep center of mass do not move
}
  



void Structureless_SystemMove()
{
	int n;

	double xcenter=0.0,ycenter=0.0;//drop moving in x and y direction to maintain the mass center of half droplets just moving in z-direction
	
	for(n=0;n<HalfNmol;n++)
	{  xcenter=xcenter+halfmol[n].r.x;
	   ycenter=ycenter+halfmol[n].r.y;
	
	}

	for (n=0;n<HalfNmol;n++)
	{  
	
		if(StepHeight==0) {halfmol[n].r.z=halfmol[n].r.z+InitialMove;}
		else { halfmol[n].r.z=halfmol[n].r.z+StepHeight+InitialMove;}
		halfmol[n].r.x=halfmol[n].r.x-xcenter/HalfNmol;//No need to do this??
		halfmol[n].r.y=halfmol[n].r.y-ycenter/HalfNmol;

	}


}

void Structureless_CenterMove()
{
	int n;

	double xcenter=0.0,ycenter=0.0;//drop moving in x and y direction to maintain the mass center of half droplets just moving in z-direction
	
	for(n=0;n<HalfNmol;n++)
	{  xcenter=xcenter+halfmol[n].r.x;
	   ycenter=ycenter+halfmol[n].r.y;
	
	}

	for (n=0;n<HalfNmol;n++)
	{  
	
	//	halfmol[n].r.z=halfmol[n].r.z+NSTEPS*stepHeight+InitialMove;
		halfmol[n].r.x=halfmol[n].r.x-xcenter/HalfNmol;
		halfmol[n].r.y=halfmol[n].r.y-ycenter/HalfNmol;

	}


}


void PrintEquilCoordVelAcc()
{
    int n;
	ofstream outfile_droplet("structurelss_halfdroplet.txt",ios::out);

    cout<<"here in print"<<endl;
//	outfile_droplet<<"Number of atoms="<<setw(14)<<HalfNmol<<endl;
	
	outfile_droplet<<"r.x"<<setw(14)<<"r.y"<<setw(14)<<"r.z"<<endl;
	
	for (n=0;n<HalfNmol;n++)
	{
		outfile_droplet<<halfmol[n].r.x<<setw(14)<<halfmol[n].r.y<<setw(14)<<halfmol[n].r.z<<endl;
	  
	 }


    outfile_droplet<<endl<<"rv.x"<<setw(14)<<"rv.y"<<setw(14)<<"rv.z"<<endl;
	for (n=0;n<HalfNmol;n++)
	{
		outfile_droplet<<halfmol[n].rv.x<<setw(14)<<halfmol[n].rv.y<<setw(14)<<halfmol[n].rv.z<<endl;
	  
	 }

    outfile_droplet<<endl<<"ra.x"<<setw(14)<<"ra.y"<<setw(14)<<"ra.z"<<endl;
	for (n=0;n<HalfNmol;n++)
	{
		outfile_droplet<<halfmol[n].ra.x<<setw(14)<<halfmol[n].ra.y<<setw(14)<<halfmol[n].ra.z<<endl;
	  
	 }

  
	outfile_droplet.close();


}




void PrintStructurelessCoordEvalution()
{
    int n;

	ofstream outfile_dropevolution("Dropinitcoord1.txt",ios::app);
	ofstream outfile_velocity("Dropinitvel1.txt",ios::app);
	ofstream outfile_acceleration("Dropinitacc1.txt",ios::app);

	for (n=0;n<HalfNmol;n++)
	{
		outfile_dropevolution<<halfmol[n].r.x<<setw(14)<<halfmol[n].r.y<<setw(14)<<halfmol[n].r.z<<endl;
		outfile_velocity<<halfmol[n].rv.x<<setw(14)<<halfmol[n].rv.y<<setw(14)<<halfmol[n].rv.z<<endl;
		outfile_acceleration<<halfmol[n].ra.x<<setw(14)<<halfmol[n].ra.y<<setw(14)<<halfmol[n].ra.z<<endl;
	 }
	outfile_dropevolution.close();
	outfile_velocity.close();
	outfile_acceleration.close();


}


void PrintStructurelessCoordEvalution2()
{
    int n;


	ofstream outfile_dropevolution2("Dropinitcoord2.txt",ios::app);	
	ofstream outfile_velocity2("Dropinitvel2.txt",ios::app);	
	ofstream outfile_acceleration2("Dropinitacc2.txt",ios::app);
	


	for (n=0;n<HalfNmol;n++)
	{
		outfile_dropevolution2<<halfmol[n].r.x<<setw(14)<<halfmol[n].r.y<<setw(14)<<halfmol[n].r.z<<endl;
		outfile_velocity2<<halfmol[n].rv.x<<setw(14)<<halfmol[n].rv.y<<setw(14)<<halfmol[n].rv.z<<endl;
		outfile_acceleration2<<halfmol[n].ra.x<<setw(14)<<halfmol[n].ra.y<<setw(14)<<halfmol[n].ra.z<<endl;
	 }
	outfile_dropevolution2.close();
	outfile_velocity2.close();
	outfile_acceleration2.close();


}



void PrintStructurelessCoordEvalution3()
{
    int n;


	ofstream outfile_dropevolution3("Dropinitcoord3.txt",ios::app);	
	ofstream outfile_velocity3("Dropinitvel3.txt",ios::app);	
	ofstream outfile_acceleration3("Dropinitacc3.txt",ios::app);
	


	for (n=0;n<HalfNmol;n++)
	{
		outfile_dropevolution3<<halfmol[n].r.x<<setw(14)<<halfmol[n].r.y<<setw(14)<<halfmol[n].r.z<<endl;
		outfile_velocity3<<halfmol[n].rv.x<<setw(14)<<halfmol[n].rv.y<<setw(14)<<halfmol[n].rv.z<<endl;
		outfile_acceleration3<<halfmol[n].ra.x<<setw(14)<<halfmol[n].ra.y<<setw(14)<<halfmol[n].ra.z<<endl;
	 }
	outfile_dropevolution3.close();
	outfile_velocity3.close();
	outfile_acceleration3.close();


}


/////////

void PrintStructurelessCoordEvalution4()
{
    int n;


	ofstream outfile_dropevolution4("Dropinitcoord4.txt",ios::app);	
	ofstream outfile_velocity4("Dropinitvel4.txt",ios::app);	
	ofstream outfile_acceleration4("Dropinitacc4.txt",ios::app);
	


	for (n=0;n<HalfNmol;n++)
	{
		outfile_dropevolution4<<halfmol[n].r.x<<setw(14)<<halfmol[n].r.y<<setw(14)<<halfmol[n].r.z<<endl;
		outfile_velocity4<<halfmol[n].rv.x<<setw(14)<<halfmol[n].rv.y<<setw(14)<<halfmol[n].rv.z<<endl;
		outfile_acceleration4<<halfmol[n].ra.x<<setw(14)<<halfmol[n].ra.y<<setw(14)<<halfmol[n].ra.z<<endl;
	 }
	outfile_dropevolution4.close();
	outfile_velocity4.close();
	outfile_acceleration4.close();


}
//////////




/////////

void PrintStructurelessCoordEvalution5()
{
    int n;


	ofstream outfile_dropevolution5("Dropinitcoord5.txt",ios::app);	
	ofstream outfile_velocity5("Dropinitvel5.txt",ios::app);	
	ofstream outfile_acceleration5("Dropinitacc5.txt",ios::app);
	


	for (n=0;n<HalfNmol;n++)
	{
		outfile_dropevolution5<<halfmol[n].r.x<<setw(14)<<halfmol[n].r.y<<setw(14)<<halfmol[n].r.z<<endl;
		outfile_velocity5<<halfmol[n].rv.x<<setw(14)<<halfmol[n].rv.y<<setw(14)<<halfmol[n].rv.z<<endl;
		outfile_acceleration5<<halfmol[n].ra.x<<setw(14)<<halfmol[n].ra.y<<setw(14)<<halfmol[n].ra.z<<endl;
	 }
	outfile_dropevolution5.close();
	outfile_velocity5.close();
	outfile_acceleration5.close();


}
//////////


/////////

void PrintStructurelessCoordEvalution6()
{
    int n;


	ofstream outfile_dropevolution6("Dropinitcoord6.txt",ios::app);	
	ofstream outfile_velocity6("Dropinitvel6.txt",ios::app);	
	ofstream outfile_acceleration6("Dropinitacc6.txt",ios::app);
	


	for (n=0;n<HalfNmol;n++)
	{
		outfile_dropevolution6<<halfmol[n].r.x<<setw(14)<<halfmol[n].r.y<<setw(14)<<halfmol[n].r.z<<endl;
		outfile_velocity6<<halfmol[n].rv.x<<setw(14)<<halfmol[n].rv.y<<setw(14)<<halfmol[n].rv.z<<endl;
		outfile_acceleration6<<halfmol[n].ra.x<<setw(14)<<halfmol[n].ra.y<<setw(14)<<halfmol[n].ra.z<<endl;
	 }
	outfile_dropevolution6.close();
	outfile_velocity6.close();
	outfile_acceleration6.close();


}
//////////


//************************End of structureless simulation******************/




void PredictorStep ()
{
  double cr[] = {19.,-10.,3.}, cv[] = {27.,-22.,7.}, div = 24., wr, wv;
  int n;

  wr = Sqr (deltaT) / div;
  wv = deltaT / div;
  for (n=0;n<HalfNmol;n++){
    halfmol[n].ro = halfmol[n].r;
    halfmol[n].rvo = halfmol[n].rv;
    PR (x);
    PRV (x);
    PR (y);
    PRV (y);
    PR (z);
    PRV (z);
    halfmol[n].ra2 = halfmol[n].ra1;
    halfmol[n].ra1 = halfmol[n].ra;
  }
}


void CorrectorStep ()
{
  double cr[] = {3.,10.,-1.}, cv[] = {7.,6.,-1.}, div = 24., wr, wv;
  int n;

  wr = Sqr (deltaT) / div;
  wv = deltaT / div;
  for (n=0;n<HalfNmol;n++) {
    CR (x);
    CRV (x);
    CR (y);
    CRV (y);
    CR (z);
    CRV (z);
  }
}



 
void HalfdropInitVels ()
{
  int n;
  vvSum=0;


   
 for (n=0;n<HalfNmol;n++)
 {
	 double sigma; 
	 sigma=sqrt(temperature);
	 VRand (&halfmol[n].rv,sigma);
	VScale (halfmol[n].rv, velMag);
	
    VVAdd (vSum, halfmol[n].rv);
	

	
  }
 
 for (n=0;n<HalfNmol;n++) {
	 VVSAdd (halfmol[n].rv, - 1. / HalfNmol, vSum);
     vvSum += VLenSq (halfmol[n].rv);}
 double temp1=(velMag*velMag*HalfNmol)/(3*(HalfNmol-1));
 cout<<temp1<<setw(14)<<vvSum/(3*(HalfNmol-1))<<endl;
 VZero (vSum);
 for (n=0;n<HalfNmol;n++) { VVAdd (vSum, halfmol[n].rv);}


 
}




/*
void HalfdropInitVels() 
{
	int n; 
	vvSum=0;
	VZero (vSum);
	for (n=0;n<HalfNmol;n++) {
		halfmol[n].rv.x=gasdev();
		halfmol[n].rv.y=gasdev();
		halfmol[n].rv.z=gasdev();

	}
	for (n=0;n<HalfNmol;n++) 
	{
		vvSum += VLenSq (halfmol[n].rv);
		VVAdd (vSum, halfmol[n].rv);
	}
	for (n=0;n<HalfNmol;n++) { VVSAdd (halfmol[n].rv, - 1. / HalfNmol, vSum);}
	double lambda = sqrt( 3 * (HalfNmol-1) * temperature / vvSum );
	cout<<"lambda" << lambda<<endl;
	for (n=0;n<HalfNmol;n++) {VScale (halfmol[n].rv, lambda);}

}


	*/


/*


void HalfdropInitVels ()
{
int n=0;
	ifstream inVel("Mrand.txt",ios::in|ios::binary);
	



	if(!inVel)
	{
		cout<<"could not open file"<<endl;
	}

	for (n=0;n<HalfNmol;n++)
	{
		

		inVel>>halfmol[n].rv.x;
		inVel>>halfmol[n].rv.y;
		inVel>>halfmol[n].rv.z;
		VScale (halfmol[n].rv, velMag);
	    VVAdd (vSum, halfmol[n].rv);


		
	}
	for (n=0;n<HalfNmol;n++) { VVSAdd (halfmol[n].rv, - 1. / HalfNmol, vSum);}
	inVel.close();
	PrinthalfdropletCoordInit();
	

}
*/

void MomentumConservation()
{
	int n;

	VZero (vSum);
	//VZero(fSum);
	for (n=0;n<HalfNmol;n++)
	{
		VVAdd (vSum, halfmol[n].rv);
		//VVAdd (fSum, halfmol[n].ra);
	}
	for (n=0;n<HalfNmol;n++) 
	{ 
		VVSAdd (halfmol[n].rv, - 1. / HalfNmol, vSum);
		//VVSAdd (halfmol[n].ra, - 1. / HalfNmol, fSum);

    }
}



void HalfdropInitAccels ()
{
  int n;

  for (n=0;n<HalfNmol;n++)
  { 
	  VZero (halfmol[n].ra);
	  VZero (halfmol[n].ra1);
	  VZero (halfmol[n].ra2);
   }

}



/*
void HalfdropInitAccels ()
{
  int n;
  double vvaSum=0;
    VecR vaSum;
   VZero (vaSum);
 for (n=0;n<HalfNmol;n++)
 {
	 //double sigma; 
	 //sigma=sqrt(temperature);
	 VRand (&halfmol[n].ra,1);
	VScale (halfmol[n].ra, -1);
	
    VVAdd (vaSum, halfmol[n].ra);
	

	
  }

}
*/


void PrintStructurelessForceEvalution(void)
{

	int n;
	ofstream outfile_dropforce("force.txt",ios::app);
	outfile_dropforce<<"time="<<setw(14)<<stepCount*deltaT<<endl;
	outfile_dropforce<<"r.z"<<setw(14)<<"ra.z"<<endl;
	for (n=0;n<HalfNmol;n++)
	{
		outfile_dropforce<<halfmol[n].r.z<<setw(14)<<halfmol[n].ra.z<<endl;
	}
  
	outfile_dropforce.close();
}



void timecounting( int tused)
{
	
	ofstream time_used("time_used.txt",ios::out);
	time_used<<"time="<<setw(14)<<tused<<endl;
	time_used.close();
}


double ContactAngleCalculation()
{
	int n,m=0;
	double ca=0;
	double R0; //initial drop radius
	double ycenter=0.0;
	int UppNmol;


	/*step 1: to cut droplet by keeping all the atom above the steps ( by requiring all StepHeight*d0<=z<=2R0)  and project all the atoms onto y-z plane)***/

	R0=initUcell.y;

	for(n=0;n<HalfNmol;n++)
	{ 
	   fitmol[n].r.x=halfmol[n].r.x;
	   fitmol[n].r.y=halfmol[n].r.y;
	   fitmol[n].r.z=halfmol[n].r.z;
	   ycenter=ycenter+fitmol[n].r.y; //center of y 
	  	
	}


	for (n=0;n<HalfNmol;n++)
	{  
	    	
		fitmol[n].r.y=fitmol[n].r.y-ycenter/HalfNmol; //moving
		
	}
	
	
	for(n=0;n<HalfNmol;n++)
	{
		if(stepCount==0) //initial step
		{
			if(fitmol[n].r.z<=2*R0)
			{
		
				
				Uppmol[m].r.y=fitmol[n].r.y;
				Uppmol[m].r.z=fitmol[n].r.z;
				m=m+1;
			}
			UppNmol=m;  //the number of atom of upper part of droplet

		}
		
		else
		{

			if(fitmol[n].r.z>=StepHeight*d0 && fitmol[n].r.z<=2*R0)
			{
		
				
				Uppmol[m].r.y=fitmol[n].r.y;
				Uppmol[m].r.z=fitmol[n].r.z;
				m=m+1;
			}
			UppNmol=m;  //the number of atom of upper part of droplet
		
		}
	
	}

	/*cout<<"UppNmol="<<UppNmol;
	for(n=0;n<30;n++)
	{
		cout<<Uppmol[n].r.y<<setw(14)<<Uppmol[n].r.z<<endl;
	}
*/

       
	//step 2: to determine the rectangular size to contain all the atoms
	double ztop=Uppmol[0].r.z, zbottom=Uppmol[0].r.z, ytop=Uppmol[0].r.y, ybottom=Uppmol[0].r.y;
	for(n=1;n<UppNmol;n++)
	{
		if(Uppmol[n].r.z<=zbottom) {zbottom=Uppmol[n].r.z;}
		if(Uppmol[n].r.z>=ztop) {ztop=Uppmol[n].r.z;}
		if(Uppmol[n].r.y<=ybottom) {ybottom=Uppmol[n].r.y;}
		if(Uppmol[n].r.y>=ytop) {ytop=Uppmol[n].r.y;}
	}


//	cout<<"ztop="<<ztop<<setw(14)<<"zbottom="<<zbottom<<setw(14)<<"ytop="<<ytop<<setw(14)<<"ybottom="<<ybottom<<endl;


   //////step 3: to determine the dropheight by N(z)/N(D)=0.9, from the bottom
	int Zlayers=int(ztop-zbottom), NFitatoms;
	double DH;
	//Fitatoms=Uppmol;
	for(n=0;n<Zlayers;n++)
	{
		NFitatoms=0;
		for(m=0;m<UppNmol;m++)
		{
			if(Uppmol[m].r.z<=(zbottom+(n+1)*(ztop-zbottom)/Zlayers))
			{
		
					//Fitatoms[NFitatoms].r.x=Uppmol[m].r.x;   
				    Fitatoms[NFitatoms].r.y=Uppmol[m].r.y;   
					Fitatoms[NFitatoms].r.z=Uppmol[m].r.z;  
					NFitatoms=NFitatoms+1;
					
			}
							
		}

	//	cout<<"n="<<n<<setw(14)<<setw(14)<<"height="<<zbottom+(n+1)*(ztop-zbottom)/Zlayers<<setw(14)<<NFitatoms<<endl;
		//for(m=0;m<NFitatoms;m++)
		//{
		//	cout<<Fitatoms[m].r.y<<setw(14)<<Fitatoms[m].r.z<<endl;
	//	}

		if(double(NFitatoms)/double(UppNmol)>=0.9)
		{
			DH=(n+1)*(ztop-zbottom)/Zlayers;
			break;
		}

	}
   // cout<<DH<<setw(14)<<Zlayers<<endl<<endl;

	
	
	/*********step 4: to determine the centeral atoms of the droplet by setting N(y)/N(Fitatomms)=0.9*/
	int NCentralDroplet;
	//double InterfaceData[100][2];
	double CenterFitatoms=0.0;
	//CentralDroplet=Fitatoms;
	for(m=0;m<NFitatoms;m++)
	{
			CenterFitatoms=CenterFitatoms+Fitatoms[m].r.y;
	}
	CenterFitatoms=CenterFitatoms/NFitatoms;
	//cout<<CenterFitatoms<<endl;
		
	for(n=0;n<1000;n++)
	{
		NCentralDroplet=0;
		for(m=0;m<NFitatoms;m++)
		{
			if(fabs(Fitatoms[m].r.y-CenterFitatoms)<=0.05+n*0.05)
			{
		
				
					CentralDroplet[NCentralDroplet].r.y=Fitatoms[m].r.y;
					CentralDroplet[NCentralDroplet].r.z=Fitatoms[m].r.z;
					NCentralDroplet=NCentralDroplet+1;
			}
			
		}
		if(double(NCentralDroplet)/double(NFitatoms)>=0.95)
		{
			break;
		}

	}

//	for(n=0;n<NCentralDroplet;n++)
//	{ 
//	    if(n<30) {cout<<CentralDroplet[n].r.y<<endl;}
	  	
//	}


	//cout<<"DH="<<double(NCentralDroplet)/double(NFitatoms)<<setw(20)<<"NCentralDroplet="<<double(NCentralDroplet)<<setw(14)<<double(NFitatoms)<<endl;

	/********************step 5:determine the discrete boundary of layered central droplet*/
	ycenter=0.0;
	for(n=0;n<NCentralDroplet;n++)
	{ 
	   ycenter=ycenter+CentralDroplet[n].r.y; //center of y 
	  	
	}
	ycenter=ycenter/NCentralDroplet;
//	cout<<"ycenter="<<ycenter<<endl;
	for(n=0;n<NCentralDroplet;n++)
	{ 
	   CentralDroplet[n].r.y=CentralDroplet[n].r.y-ycenter; //center of y 
	  // if(n<30) {cout<<CentralDroplet[n].r.y<<endl;}
	  	
	}
	

	//int NFitPoints=100; //using 100 points to fit interface;
	int q,s,t,mcenter=0;
	double extremity=0.0;
	//LayerAtoms=CentralDroplet;

	

	if(DH<2.0) { ca=0.0;goto end;} //if the height is very small, we think it is complete wetting
	//cout<<"DH="<<DH<<setw(28)<<"NCentralDroplet="<<NCentralDroplet<<endl;

	for(n=0;n<NFitPoints;n++)//layers
	{
		m=0;
		for(q=0;q<NCentralDroplet;q++)
		{ 
			if(CentralDroplet[q].r.z>=zbottom+n*DH/NFitPoints && CentralDroplet[q].r.z<zbottom+(n+1)*DH/NFitPoints)
			{
				LayerAtoms[m].r.y=CentralDroplet[q].r.y; //center of y 
				LayerAtoms[m].r.z=CentralDroplet[q].r.z;
				m=m+1;  //the number of atom of nth layer
			}
						
		}
	
	//	cout<<"n="<<n<<setw(14)<<endl;
	//	cout<<"m="<<m<<setw(14)<<endl;

		if(m<=0) {cout<<"wrong in determine the points on the interface"<<endl;exit(0);}  //
		ycenter=0;
		for(t=0;t<m;t++)
		{

			ycenter=ycenter+LayerAtoms[t].r.y; 
		/*	if(n==0)
			{
				cout<<LayerAtoms[t].r.y<<endl;
			}*/
		}
		
	    ycenter=ycenter/m; //to decide the center of the nth layer
		//cout<<n<<" th layer"<<setw(14)<<ycenter<<endl; good
		/*for(t=0;t<m;t++)
		{
			LayerAtoms[t].r.y=LayerAtoms[t].r.y-ycenter; //to shift the atoms to the center of the nth layer

			

		}*/
		
		
		
		for(t=0;t<2000;t++)
		{
			mcenter=0; //to count the number of atom within center of nth layer by mcenter/m <=0.95
			for(s=0;s<m;s++)
			{
				if(fabs(LayerAtoms[s].r.y-ycenter)<=0.05+0.05*t)
				{
					mcenter=1+mcenter;
					
				}
			}
			if((double(mcenter)/double(m))>=0.95) {extremity=0.05+t*0.05;break;}

		}
	//	cout<<"m="<<m<<setw(14)<<"mcenter="<<mcenter<<endl;

		InterfaceData[n][0]=extremity;
		InterfaceData[n][1]=zbottom+(n+1)*DH/NFitPoints;
	}

	/*for(n=0;n<NFitPoints;n++)
	{
		cout<<"n=:"<<n<<setw(14)<<InterfaceData[n][0]<<setw(14)<<InterfaceData[n][1]<<endl;
	}*/


   /*********step 6 to determine the the center, and the radius of fitting circle **************/
	circlefit(); //return xc,yc, and rc
//	cout<<"xc="<<circle[0]<<setw(14)<<"yc="<<circle[1]<<setw(14)<<"rc="<<circle[2]<<endl; 

   /*to determine the contact angle*********/

	double dydx, xb,yb;
	yb=zbottom;
	xb=sqrtl(circle[2]*circle[2]-(yb-circle[1])*(yb-circle[1]))+circle[0];
	//cout<<"xb="<<xb<<endl;
//	cout<<"yb="<<yb<<endl;

	dydx=-(xb-circle[0])/(yb-circle[1]);
	//cout<<"dydx="<<dydx<<endl;
	if(dydx>0){ ca=180-atanl(dydx)*180/M_PI;}
	else{ ca=-atanl(dydx)*180/M_PI;}

    end: 
	return (ca);

}









void circlefit()
{
	//circle fit by Pratt V.Pratt, "Direct least-squares fitting of algebraic surfaces", computer Graphics, v 21, p145-152(1987) 
	
	
	///////////////////to store all fitting data//////////////
	int m,n,k,s;
	double XY[2*NFitPoints][2]; 
	for(m=0;m<2*NFitPoints;m++)
	{
		if(m<NFitPoints)
		{
			XY[m][0]=InterfaceData[m][0];
			XY[m][1]=InterfaceData[m][1];
		}
		if(m>=NFitPoints)
		{
			XY[m][0]=-InterfaceData[m-NFitPoints][0];
			XY[m][1]=InterfaceData[m-NFitPoints][1];
		}
	}

   ////////to form matrix A and B such that A x=B, where X=(xc,yc,rc)//////////////
	double A[2*NFitPoints][3],B[2*NFitPoints];
	double d;
	double xc=0,yc=0,rc=initUcell.x;  //initial guess
	int NIter=100; //Maximal iteration for Gauss-Newton method

	for(s=0;s<NIter;s++)
	{
		for(m=0;m<2*NFitPoints;m++)
		{
		    d=sqrtl((XY[m][0]-xc)*(XY[m][0]-xc)+(XY[m][1]-yc)*(XY[m][1]-yc));
			A[m][0]=(xc-XY[m][0])/d;
			A[m][1]=(yc-XY[m][1])/d;
			A[m][2]=-1;
			B[m]=-(d-rc);
		
		}

	
		////to get ATATranspose[A].[A] and transpose[A].B (by least-square method
	//	double ATA[3][3],ATB[3];
		for(m=0;m<3;m++)
		{
			for(n=0;n<3;n++)
			{
				ATA[m][n]=0.0;
				for(k=0;k<2*NFitPoints;k++)
				{
					ATA[m][n]=ATA[m][n]+A[k][m]*A[k][n];
				}
			}
			ATB[m]=0.0;
			for(k=0;k<2*NFitPoints;k++)
			{
				ATB[m]=ATB[m]+A[k][m]*B[k];
			}
		}

    
    	
		///to solve ATATranspose[A].[A]*x=transpose[A].B
		m=aldle();
		if(sqrtl(dcircle[0]*dcircle[0]+dcircle[1]*dcircle[1]+dcircle[2]*dcircle[2])<epsilon)
		{
			circle[0]=dcircle[0]+xc;
			circle[1]=dcircle[1]+yc;
			circle[2]=dcircle[2]+rc;
			break;
		}
		else
		{
			//cout<<"s="<<s<<endl;
			//cout<<"error=="<<sqrtl(dcircle[0]*dcircle[0]+dcircle[1]*dcircle[1]+dcircle[2]*dcircle[2])<<endl;
			if(s==NIter-1){cout<<"No solution to search a circle by Gauss-Newton method"<<endl;break;}
		//	cout<<"xc="<<xc<<setw(14)<<"yc="<<yc<<setw(14)<<"rc="<<rc<<endl;
			xc=dcircle[0]+xc;
			yc=dcircle[1]+yc;
			rc=dcircle[2]+rc;
		}


	}
}

void PrintInterface()
{

	int m;
	ofstream outfile_interface("interface.txt",ios::app);
	outfile_interface<<"time="<<setw(14)<<stepCount*deltaT<<endl;
	outfile_interface<<"y"<<setw(14)<<"z"<<endl;
	for(m=0;m<NFitPoints;m++)
	{
		outfile_interface<<InterfaceData[m][0]<<setw(14)<<InterfaceData[m][1]<<endl;
	 
	 }
  
	outfile_interface.close();


}

void  PrintContactAngle(double ca)
{
	
	
	ofstream outfile_contactangle("ContactAngle.txt",ios::app);
	outfile_contactangle<<stepCount*deltaT<<setw(14)<<ca<<endl;
	outfile_contactangle.close();
}


//The goal of this rountine is to solve Ax=b at A=Transpose[T]. Input informtion is need of A and b ,and solution are stored in b
int aldle()
{
	double A[3][3],B[3];
	double L[3][3];
	double D[3][3];
	double y[3],x[3];
	int m,n,M=3,N=3;
	//to decompose A=L*D*Transpose[L]
	
	for(m=0;m<M;m++)
	{
		for(n=0;n<N;n++)
		{
			A[m][n]=ATA[m][n];
			if(m==n){ L[m][n]=1;D[m][n]=1;} //initialize L and D matrice
			else {L[m][n]=0; D[m][n]=0;}
		}
		B[m]=ATB[m];
	}

	
	//to iterate to get L and D (*checked*)
	int k;
	double sum1,sum2;
	D[0][0]=A[0][0];
	for(m=1;m<M;m++)
	{
		
		
		for(n=0;n<m;n++)
		{
			sum2=0;
			for(k=0;k<=n-1;k++)
			{
				sum2=sum2+L[m][k]*L[n][k]*D[k][k];
			}
		
			L[m][n]=(A[m][n]-sum2)/D[n][n];
		}
		
		sum1=0;
		for(k=0;k<=m-1;k++)
		{
			sum1=sum1+L[m][k]*L[m][k]*D[k][k];
		}
		D[m][m]=A[m][m]-sum1;
		
		
	}

   // to get solution of x 
	y[0]=B[0];
	for(m=1;m<M;m++)
	{
		sum1=0;
		for(k=0;k<=m-1;k++)
		{
			sum1=sum1+L[m][k]*y[k];
		}
		y[m]=B[m]-sum1;
	}
	

	x[M-1]=y[M-1]/D[M-1][M-1];

	for(m=M-1;m>=0;m--)
	{
		sum2=0;
		for(k=m+1;k<M;k++)
		{
			sum2=sum2+D[m][m]*L[k][m]*x[k];
		}
		x[m]=(y[m]-sum2)/D[m][m];
	}
	for(m=0;m<M;m++)
	{
		dcircle[m]=x[m];//to return dxc,dyc, and drc and stored in dcircle[3]
	}

	

	return (1);
}



void PrinthalfdropletCoordInit()
{

	int n;
	
	ofstream outfile_halfdroplet("halfdropletinit.txt",ios::out);
	

	outfile_halfdroplet<<"HalfNmol="<<HalfNmol<<endl;
	//outfile_halfdroplet<<"r.x"<<setw(14)<<"r.y"<<setw(14)<<"r.z"<<endl;
	//for (n=0;n<HalfNmol;n++)
	//{
	//	outfile_halfdroplet<<halfmol[n].r.x<<setw(14)<<halfmol[n].r.y<<setw(14)<<halfmol[n].r.z<<endl;
	// }
	outfile_halfdroplet<<"r.vx"<<setw(14)<<"r.vy"<<setw(14)<<"r.z"<<endl;
	for (n=0;n<HalfNmol;n++)
	{
		outfile_halfdroplet<<halfmol[n].rv.x<<setw(14)<<halfmol[n].rv.y<<setw(14)<<halfmol[n].rv.z<<endl;

	 }
//	outfile_halfdroplet<<"r.ax"<<setw(14)<<"r.ay"<<setw(14)<<"r.z"<<endl;
	//for (n=0;n<HalfNmol;n++)
//	{
	//	outfile_halfdroplet<<halfmol[n].ra.x<<setw(14)<<halfmol[n].ra.y<<setw(14)<<halfmol[n].ra.z<<endl;

	// }




/*	outfile_halfdroplet<<"rv.x"<<setw(14)<<"rv.y"<<setw(14)<<"rv.z"<<endl;
	for (n=0;n<HalfNmol;n++)
	{
		outfile_halfdroplet<<halfmol[n].rv.x<<setw(14)<<halfmol[n].rv.y<<setw(14)<<halfmol[n].rv.z<<endl;
	  
	 }

	outfile_halfdroplet<<"ra.x"<<setw(14)<<"ra.y"<<setw(14)<<"ra.z"<<endl;
	for (n=0;n<HalfNmol;n++)
	{
		outfile_halfdroplet<<halfmol[n].ra.x<<setw(14)<<halfmol[n].ra.y<<setw(14)<<halfmol[n].ra.z<<endl;
	  
	 }*/
    
	outfile_halfdroplet.close();
}


void forceztotal(void)
{

	int n;
	ofstream outfile_dropforce("force.txt",ios::app);
	double fx=0.0,fy=0.0,fz=0.0;
	//outfile_dropforce<<"time="<<setw(14)<<stepCount*deltaT<<endl;
	for (n=0;n<HalfNmol;n++)
	{
		fx=fx+halfmol[n].ra.x;
		fy=fy+halfmol[n].ra.y;
		fz=fz+halfmol[n].ra.z;
	}

	outfile_dropforce<<stepCount*deltaT<<setw(14)<<fx/HalfNmol<<setw(14)<<fy/HalfNmol<<setw(14)<<fz/HalfNmol<<endl;
  
	outfile_dropforce.close();
}


void PrintVelocityDistribution()
{
	int n,m;
	double vxmin=0.0,vxmax=0.0,vymin=0.0,vymax=0.0, vzmin=0.0,vzmax=0.0;
	double deltax,deltay,deltaz;
	double nx,ny,nz;
	double nxx=0,nyy=0,nzz=0;

    for (n=0;n<HalfNmol;n++)
	{
		if(halfmol[n].rv.x<vxmin){vxmin=halfmol[n].rv.x;}
		if(halfmol[n].rv.x>=vxmax){vxmax=halfmol[n].rv.x;}

		if(halfmol[n].rv.y<vymin){vymin=halfmol[n].rv.y;}
		if(halfmol[n].rv.y>=vymax){vymax=halfmol[n].rv.y;}

		if(halfmol[n].rv.z<vzmin){vzmin=halfmol[n].rv.z;}
		if(halfmol[n].rv.z>=vzmax){vzmax=halfmol[n].rv.z;}

	 }
	deltax=(vxmax-vxmin)/100;
	deltay=(vymax-vymin)/100;
	deltaz=(vzmax-vzmin)/100; //output 1000 points for velocity distribution


	ofstream outfile_vd("velocity distribution.txt",ios::app);
	outfile_vd<<"timenow="<<timeNow<<endl;
	outfile_vd<<"px"<<setw(14)<<"rv.x"<<setw(14)<<"py"<<setw(14)<<"rv.y"<<setw(14)<<"pz"<<setw(14)<<"rv.z"<<endl;
	for(m=0;m<101;m++)
	{
		nx=0.;
		ny=0.;
		nz=0.;
		for(n=0;n<HalfNmol;n++)
		{
			
			if((m*deltaz+vzmin<=halfmol[n].rv.z) && (halfmol[n].rv.z <(m+1)*deltaz+vzmin)) {nz=nz+1.;}
			if((m*deltay+vymin<=halfmol[n].rv.y) && (halfmol[n].rv.y <(m+1)*deltay+vymin)) {ny=ny+1.;}
			if((m*deltax+vxmin<=halfmol[n].rv.x) && (halfmol[n].rv.x <(m+1)*deltax+vxmin)) {nx=nx+1.;}
			
			
		}
	
		outfile_vd<<nx/HalfNmol<<setw(14)<<(m+0.5)*deltax+vxmin<<setw(14)        \
			      <<ny/HalfNmol<<setw(14)<<(m+0.5)*deltay+vymin<<setw(14)                      \
			      <<nz/HalfNmol<<setw(14)<<(m+0.5)*deltaz+vzmin<<endl;
	}

	outfile_vd.close();




 }

double AForce(double yl, double ys, double zl, double zs)
{
	double fcVal, dur, dua;
	if(fabs(yl)<tol)
	{
		dur=AFUrep(yl,zl)+FUrep(ys,zs)-FUrep(ys,zl)-AFUrep(yl,zs);
		dua=AFUatt(yl,zl)+FUatt(ys,zs)-FUatt(ys,zl)-AFUatt(yl,zs);
		fcVal=4.0*Egs*rous*(Lgs12*dur-Lgs6*dua);
	}

	if(fabs(ys)<tol)
	{
		dur=FUrep(yl,zl)+AFUrep(ys,zs)-AFUrep(ys,zl)-FUrep(yl,zs);
		dua=FUatt(yl,zl)+AFUatt(ys,zs)-AFUatt(ys,zl)-FUatt(yl,zs);
		fcVal=4.0*Egs*rous*(Lgs12*dur-Lgs6*dua);
	}
	if(fabs(zl)<tol ||fabs(zs)<tol)
	{

		dur=FUrep(yl,zl)+FUrep(ys,zs)-FUrep(ys,zl)-FUrep(yl,zs);
		dua=FUatt(yl,zl)+FUatt(ys,zs)-FUatt(ys,zl)-FUatt(yl,zs);
		fcVal=4.0*Egs*rous*(Lgs12*dur-Lgs6*dua);

	}
	return(fcVal);

}

double FUatt(double y,double z)//////
{
	double fcVal;
	fcVal=M_PI*z*(3*y*y+2*z*z)/(8.*powl(y,4.0)*powl(y*y+z*z,1.5));
	return(fcVal);
}
double FUrep(double y,double z)
{
	double fcVal;
	fcVal=M_PI*z*(315*powl(y,8.0)+840*powl(y,6)*z*z+1008*powl(y*z,4.0)+576*y*y*powl(z,6.0)+128*powl(z,8.0))/(1280.*powl(y,10.0)*powl(y*y+z*z,4.5));
	return(fcVal);
}

double AFUatt(double y, double z)
{
	double fcVal;
	fcVal=Sign(z)*M_PI*(-3/(32*powl(z,4.0))+5.*y*y/(32*powl(z,6.0))-105*powl(y,4.0)/(512*powl(z,8.0)));
	return(fcVal);
}

double AFUrep(double y, double z)
{
	double fcVal;
	fcVal=Sign(z)*M_PI*(-63/(2560.*powl(z,10.0))+231.*y*y/(2048*powl(z,12.0))-1287*powl(y,4.0)/(4096*powl(z,14.0)));
	return(fcVal);
}



