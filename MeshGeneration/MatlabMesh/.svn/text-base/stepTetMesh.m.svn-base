function [Nodes,Elements,Boundary]=stepTetMesh(level)

  %% Meshes channel into boxes (except in step), then splits boxes into 6
  %% tetrahedra following picture on Maubach's webpage.

  %% Number of triangles on 
  %% [zmin,zmax,ymin,ymax,xmin,xmax, step_y=5, step_y=6, step_z=1]
  
  bdyidx = [0,0,0,0,0,0, 0,0,0];  

  %% 0 = Dirichlet, 1 = Neumann, 2 = "do nothing"

  bdycode = [0,0,0,2,0,0, 0,0,0];   

  if level==0
    xx = linspace(0,10,3);
    yy = [0,5,6,10,20,30,40,50];
    zz = [0,1,5,10];
  elseif level==1
    xx = linspace(0,10,5);
    yy = [0,2.5,5,5.5,6,8,10,15,20,30,40,45,50];
    zz = [0,0.5,1,2,3,6,9,10];
  elseif level==2
    xx = [0,1,3,5,7,9,10];
    yy = [0,1,3,5,5.5,6,7,8,10,12,15,18,21,25,28,31,35,40,45,50];
    zz = [0,0.33,0.66,1,1.5,2,3,5,7,9,10];    
  else
    error('invalid mesh level for step problem')
  end

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
          xco = Nodes(index,1);
          yco = Nodes(index,2);
          zco = Nodes(index,3);
	  %% Get the 8 nodes corresponding to this cube
          BRF = BLF + 1;
          BLB = BLF + N1;
          BRB = BLB + 1;
          TLF = BLF + N1*N2;
          TRF = BRF + N1*N2;
          TLB = BLB + N1*N2;
          TRB = BRB + N1*N2;
          
	  %% divide cube into 6 tetrahedra
          tetrahedra(1,1:4)=[TRF,BLF,BRF,BLB];
          tetrahedra(2,1:4)=[BRB,TRF,BRF,BLB];
          tetrahedra(3,1:4)=[TRB,TRF,BRB,BLB];
          tetrahedra(4,1:4)=[TLB,TRF,TRB,BLB];
          tetrahedra(5,1:4)=[TLF,TRF,TLB,BLB];
          tetrahedra(6,1:4)=[BLF,TRF,TLF,BLB];
            
          %% the cube is now built
            
          c_of_cube = 1/2 * (Nodes(BLF,:) + Nodes(TRB,:) );
             
          if c_of_cube(2)<5.99999 && c_of_cube(2)>5.00001 && c_of_cube(3)<.999999
	    %% do nothing, you are in the step
          else             
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

            if Nodes(TLF,2) < 6.000001 && Nodes(TLF,3)< 1.000001 
	      %% front is step face y=6
	      bdyidx(8) = bdyidx(8)+1;
	      Boundary(8,bdyidx(8),1:3) = [BRF, TRF, BLF];
	      bdyidx(8) = bdyidx(8)+1;
	      Boundary(8,bdyidx(8),1:3) = [TLF, BLF, TRF];
	    end

            if Nodes(TLB,2) > 4.99999 && Nodes(TLB,3)< 1.000001 
	      %% back is step face y=5
	      bdyidx(7) = bdyidx(7)+1;
	      Boundary(7,bdyidx(7),1:3) = [TLB, TRB, BLB];
	      bdyidx(7) = bdyidx(7)+1;
	      Boundary(7,bdyidx(7),1:3) = [BRB, BLB, TRB];
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

            if c_of_cube(2)<5.99999 && c_of_cube(2)>5.00001 && Nodes(BRF,3)<1.000001 %% bottom is step face z=1.0
	      bdyidx(9) = bdyidx(9)+1;
	      Boundary(9,bdyidx(9),1:3) = [BRB, BRF, BLB];
	      bdyidx(9) = bdyidx(9)+1;
	      Boundary(9,bdyidx(9),1:3) = [BLF, BLB, BRF];
	    end

          end        
        end
      end
    end
  end



% %Plot the mesh!
% 
% First draw the elements
%Nodes1=Nodes;
%for i=1:size(Elements,1)
% for i=10:50:180
%    lcn=Elements(i,:);
%     xx = zeros(20,3);
%     for j=1:20
%         xx(j,:) = Nodes1(lcn(j),:);
%         plot3(xx(j,1),xx(j,2),xx(j,3),'r*')
%         hold on
%     end
%    x1 = Nodes1(lcn(1),:);
%    x2 = Nodes1(lcn(2),:);
%    x3 = Nodes1(lcn(3),:);
%    x4 = Nodes1(lcn(4),:);
%     plot3(x1(1),x1(2),x1(3),'r*')
%     hold on
%     plot3(x2(1),x2(2),x2(3),'r*')
%     hold on
%     plot3(x3(1),x3(2),x3(3),'r*')
%     hold on
%     plot3(x4(1),x4(2),x4(3),'r*')
%     hold on
    
%   plot3(linspace(x1(1),x2(1)),linspace(x1(2),x2(2)),linspace(x1(3),x2(3)))
%   hold on
%   plot3(linspace(x1(1),x3(1)),linspace(x1(2),x3(2)),linspace(x1(3),x3(3)))
%   hold on
%   plot3(linspace(x1(1),x4(1)),linspace(x1(2),x4(2)),linspace(x1(3),x4(3)))
%   hold on
%    plot3(linspace(x2(1),x3(1)),linspace(x2(2),x3(2)),linspace(x2(3),x3(3)))
%    hold on
%    plot3(linspace(x2(1),x4(1)),linspace(x2(2),x4(2)),linspace(x2(3),x4(3)))
%    hold on
%    plot3(linspace(x3(1),x4(1)),linspace(x3(2),x4(2)),linspace(x3(3),x4(3)))
%    hold on
% end
% zlabel('z')
% ylabel('y')
% xlabel('x')
% view(90,0)
        
  fname = sprintf("stepMesh%5.3fMesh.tplc", yy(2)-yy(1));
  fid = fopen(fname, "w");
  fprintf(fid, "TPLC\nDIMENSION 3\n");
  
  fprintf(fid, "POINTS %d\n", size(Nodes,1));
  for i=1:size(Nodes,1) 
    fprintf(fid, "%d %12.6f %12.6f %12.6f\n", i, Nodes(i,1), Nodes(i,2), Nodes(i,3));
  end

  fprintf(fid, "SEGMENTS 0\n");

  fprintf(fid, "FACES 9\n");
  for i=1:9
    fprintf(fid, "%d %d\n", 200+bdycode(i), bdyidx(i));

    for j=1:bdyidx(i) 
      fprintf(fid, "%d %d %d\n", Boundary(i,j,1), Boundary(i,j,2), Boundary(i,j,3));
    end 
  end

  fprintf(fid, "TETS %d\n", size(Elements,1));
  for i=1:size(Elements,1)
    fprintf(fid, "%d %d %d %d\n", Elements(i,1), Elements(i,2), Elements(i,3),  Elements(i,4));
  end

  fclose(fid);
