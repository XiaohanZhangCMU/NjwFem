function [Nodes,Elements,Boundary]=cubeTetMesh(level)

level = 4;

  %% Meshes [0,1] Cube into boxes then each box into 6 tets

  %% Number of triangles on [zmin,zmax,ymin,ymax,xmin,xmax]
  
  bdyidx = [0,0,0,0,0,0];  
  
  %% [-z +z -y +y -x +x]

  bdycode = [0,1,1,1,1,1];   %% 0 = Dirichlet, 1 = Neumann or "nothing"
    
  xx = (2/2^level)*[0:2^level]-1;   %% [0,1] -> [-1,1]
 
  %% XIAOHAN
 
  xx = 25*xx;
  
  %%
  yy=xx;		
  zz=xx;
  
  N1 = length(xx);
  N2 = length(yy);
  N3 = length(zz);
    
  index=0;
  for k=1:N3
    for j=1:N2
      for i=1:N1
        index=index+1;
        Nodes(index,:)=[xx(i),yy(j),zz(k)]; 
      end
    end  
  end
	  
  cubeNumber=0;
  index = 0;

  for k=1:N3
    for j=1:N2
      for i=1:N1
        index=index+1;
        if i~=N1 && j~=N2 && k~=N3
          %% BLF is the bottom left front of each cube            
          BLF = index;
          BRF = BLF + 1;
          BLB = BLF + N1;
          BRB = BLB + 1;
          TLF = BLF + N1*N2;
          TRF = BRF + N1*N2;
          TLB = BLB + N1*N2;
          TRB = BRB + N1*N2;
          
          %% divide cube into 6 tetrahedra
          
          tetrahedra(1,1:4)=[BLF,TRF,BRF,BLB];
          tetrahedra(2,1:4)=[TRF,BRB,BRF,BLB];
          tetrahedra(3,1:4)=[TRF,TRB,BRB,BLB];
          tetrahedra(4,1:4)=[TRF,TLB,TRB,BLB];
          tetrahedra(5,1:4)=[TRF,TLF,TLB,BLB];
          tetrahedra(6,1:4)=[TRF,BLF,TLF,BLB];          

         %% tetrahedra(1,1:4)=[TRF,BLF,BRF,BLB];
         %% tetrahedra(2,1:4)=[BRB,TRF,BRF,BLB];
         %% tetrahedra(3,1:4)=[TRB,TRF,BRB,BLB];
         %% tetrahedra(4,1:4)=[TLB,TRF,TRB,BLB];
         %% tetrahedra(5,1:4)=[TLF,TRF,TLB,BLB];
         %% tetrahedra(6,1:4)=[BLF,TRF,TLF,BLB];          
          
          %% the cube is now built
          
          cubeNumber=cubeNumber+1;          
          elnum=6*(cubeNumber-1)+1;
          Elements(elnum,:)=tetrahedra(1,:);
          Elements(elnum+1,:)=tetrahedra(2,:);
          Elements(elnum+2,:)=tetrahedra(3,:);
          Elements(elnum+3,:)=tetrahedra(4,:);
          Elements(elnum+4,:)=tetrahedra(5,:);
          Elements(elnum+5,:)=tetrahedra(6,:);
	  
          %% construct the boundary triangulation
	  
	  if i==1  %% Left
	    bdyidx(5) = bdyidx(5)+1;
	    Boundary(5,bdyidx(5),1:3) = [BLF, TLF, BLB];
	    bdyidx(5) = bdyidx(5)+1;
	    Boundary(5,bdyidx(5),1:3) = [TLB, BLB, TLF];
	  end
	  if i==N1-1  %% Right
	    bdyidx(6) = bdyidx(6)+1;
	    Boundary(6,bdyidx(6),1:3) = [BRF, BRB, TRF];
	    bdyidx(6) = bdyidx(6)+1;
	    Boundary(6,bdyidx(6),1:3) = [TRB, TRF, BRB];
	  end

	  if j==1  %% Front
	    bdyidx(3) = bdyidx(3)+1;
	    Boundary(3,bdyidx(3),1:3) = [BRF, TRF, BLF];
	    bdyidx(3) = bdyidx(3)+1;
	    Boundary(3,bdyidx(3),1:3) = [TLF, BLF, TRF];
	  end
	  if j==N2-1  %% Back
	    bdyidx(4) = bdyidx(4)+1;
	    Boundary(4,bdyidx(4),1:3) = [TLB, TRB, BLB];
	    bdyidx(4) = bdyidx(4)+1;
	    Boundary(4,bdyidx(4),1:3) = [BRB, BLB, TRB];
	  end

	  if k==1  %% Bottom
	    bdyidx(1) = bdyidx(1)+1;
	    Boundary(1,bdyidx(1),1:3) = [BRB, BRF, BLB];
	    bdyidx(1) = bdyidx(1)+1;
	    Boundary(1,bdyidx(1),1:3) = [BLF, BLB, BRF];
	  end
	  if k==N3-1  %% Top
	    bdyidx(2) = bdyidx(2)+1;
	    Boundary(2,bdyidx(2),1:3) =  [TRB, TLB, TRF];
	    bdyidx(2) = bdyidx(2)+1;
	    Boundary(2,bdyidx(2),1:3) =  [TLF, TRF, TLB];
	  end

        end
      end
    end
  end

%%Plot the mesh!
%%for i=1:size(Elements,1)
%%    lcn=Elements(i,:);
%%    x1 = Nodes(lcn(1),:);
%%    x2 = Nodes(lcn(2),:);
%%    x3 = Nodes(lcn(3),:);
%%    x4 = Nodes(lcn(4),:);
%%   plot3(linspace(x1(1),x2(1)),linspace(x1(2),x2(2)),linspace(x1(3),x2(3)))
%%   hold on
%%   plot3(linspace(x1(1),x3(1)),linspace(x1(2),x3(2)),linspace(x1(3),x3(3)))
%%   hold on
%%   plot3(linspace(x1(1),x4(1)),linspace(x1(2),x4(2)),linspace(x1(3),x4(3)))
%%   hold on
%%    plot3(linspace(x2(1),x3(1)),linspace(x2(2),x3(2)),linspace(x2(3),x3(3)))
%%    hold on
%%    plot3(linspace(x2(1),x4(1)),linspace(x2(2),x4(2)),linspace(x2(3),x4(3)))
%%    hold on
%%    plot3(linspace(x3(1),x4(1)),linspace(x3(2),x4(2)),linspace(x3(3),x4(3)))
%%    hold on
%%end
%%zlabel('z')
%%ylabel('y')
%%xlabel('x')
%%view(0,90)

  fname = sprintf('cubeMesh%5.3f.tplc', xx(2)-xx(1));
  fid = fopen(fname, 'w');
  fprintf(fid, 'TPLC\nDIMENSION 3\n');
  
  fprintf(fid, 'POINTS %d\n', size(Nodes,1));
  for i=1:size(Nodes,1) 
    fprintf(fid, '%d %12.6f %12.6f %12.6f\n', i, Nodes(i,1), Nodes(i,2), Nodes(i,3));
  end

  fprintf(fid, 'SEGMENTS 0\n');

  fprintf(fid, 'FACES 6\n');
  for i=1:6
    fprintf(fid, '%d %d\n', 200+bdycode(i), bdyidx(i));

    for j=1:bdyidx(i) 
      fprintf(fid, '%d %d %d\n', Boundary(i,j,1), Boundary(i,j,2), Boundary(i,j,3));
    end 
  end

  fprintf(fid, 'TETS %d\n', size(Elements,1));
  for i=1:size(Elements,1)
    fprintf(fid, '%d %d %d %d\n', Elements(i,1), Elements(i,2), Elements(i,3),  Elements(i,4));
  end

  fclose(fid);
 
