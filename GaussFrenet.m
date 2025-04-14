function [Lk_all] = GaussFrenet(filename,specialChainID,circBool,allPlotBool,targetAtom,plotBaseNum)
%copyright Simon Vecchioni, August 2023 (sv1091@nyu.edu)
%Computes various geometric and topological properties of macromolecular
     % chains in PDB files
%filename is PDB file directory 'PDBs', i.e.: "2T3_merged.pdb"
%specialChainID: which pdb chain is being analyzed? 'D', etc.
%circBool: artificially connect first and last atom? 0/1 for N/Y
%allPlotBool: show all data? 0/1 for N/Y
%targetAtom: name of atom in PDB, typically 'O3''' or 'P'
%plotBaseNum: number of nucleotide steps desired in graphs, 0 is all:
%      corresponds to "length of center strand" for some applications. For
%      example, a 14 bp triangle edge should have an asymmetric unit in the
%      pdb file with 14 bp, so set this to 14; if you want to see the whole
%      triangle/chain in the PDB file, set this to 0

[~, fname, ~] = fileparts(filename); %for saving metadata in figs
filename_load = strcat('PDBs/',filename); %change this if you change the folder name of PDBs

%make targetAtom all a cell of strings; iterate all relevant atoms
if iscell(targetAtom)

else
    targetAtom = {targetAtom};
end

currStruct = pdbread(filename_load);
model = currStruct.Model;
try
    atoms = model.Atom;
catch
    atoms = model.HeterogenAtom;
end
atomnum = size(atoms, 2);
aNames = {atoms(1:atomnum).AtomName}.';
resIDs = {atoms(1:atomnum).resSeq}.';
resNames = {atoms(1:atomnum).resName}.';
chainIDs = {atoms(1:atomnum).chainID}.';
r = zeros(3,1); %coordinates
rS = zeros(3,1); %base coordinates
P = 0; %phosphate counter
S = 0; %sugar counter
numAtomTypes = max(size(targetAtom));

%get N1/N9 for D/L-DNA; modify here if using non-standard residues
for i = 1:atomnum
    % i counts the atom; resnum counts residue
    if S < resIDs{i,1} && chainIDs{i,1}
        if chainIDs{i,1} == specialChainID
            if (strcmp('DT',resNames{i,1}) && strcmp('N1',aNames{i,1})) || ...
                    (strcmp('DC',resNames{i,1}) && strcmp('N1',aNames{i,1})) || ...
                    (strcmp('DG',resNames{i,1}) && strcmp('N9',aNames{i,1})) ||...
                    (strcmp('DA',resNames{i,1}) && strcmp('N9',aNames{i,1})) ||...
                    (strcmp('0DT',resNames{i,1}) && strcmp('N1',aNames{i,1})) ||...
                    (strcmp('0DC',resNames{i,1}) && strcmp('N1',aNames{i,1})) ||...
                    (strcmp('0DG',resNames{i,1}) && strcmp('N9',aNames{i,1})) ||...
                    (strcmp('0DA',resNames{i,1}) && strcmp('N9A',aNames{i,1}))
                S = S+1;
                rS(S,1) = atoms(i).X.';
                rS(S,2) = atoms(i).Y.';
                rS(S,3) = atoms(i).Z.';
            end
        end

    end

end


%get atoms that are relevent
for i = 1:atomnum
    for j = 1:numAtomTypes
        if strcmp(aNames{i,1},targetAtom{j,1})...
                && chainIDs{i,1} == specialChainID
            P = P+1;
            r(P,1) = atoms(i).X.';
            r(P,2) = atoms(i).Y.';
            r(P,3) = atoms(i).Z.';
        end
    end
end

if allPlotBool
    rAll = zeros(3,1); %coordinates
    PAll = 0; %phosphate counter
    for i = 1:atomnum
        for j = 1:numAtomTypes
            if strcmp(aNames{i,1},targetAtom{j,1})
                PAll = PAll+1;
                rAll(PAll,1) = atoms(i).X.';
                rAll(PAll,2) = atoms(i).Y.';
                rAll(PAll,3) = atoms(i).Z.';
            end
        end
    end
end


%% tangent vectors
tV = zeros(3,P);
for i = 1:P-1
    tV(:,i) = (r(i+1,:)-r(i,:))/euclidDist(r(i+1,:),r(i,:));
end

if circBool
    tV(:,P) = (r(1,:)-r(end,:))/euclidDist(r(1,:),r(end,:));
end

%% b vectors

bV = zeros(3,P);
for i = 2:P
    bV(:,i) = cross(tV(:,i-1)',tV(:,i)')/euclidDist(cross(tV(:,i-1)',tV(:,i)'),[0, 0 ,0]);
end

if circBool
    bV(:,1) = cross(tV(:,P)',tV(:,1)')/euclidDist(cross(tV(:,P)',tV(:,1)'),[0, 0 ,0]);
end

%% n vectors
nV = zeros(3,P);
for i = 1:P
    nV(:,i) = cross(bV(:,i),tV(:,i));
end

%% c vectors
cV = zeros(size(tV));
for i = 1:P-1
    cV(:,i) = (rS(i,:)-r(i,:))/euclidDist(rS(i,:),r(i,:));
end

if circBool
    cV(:,P) = (rS(end,:)-r(end,:))/euclidDist(rS(end,:),r(end,:));
end

%% p vectors
pV = zeros(3,P);
for i = 2:P
    pV(:,i) = cross(tV(:,i),cV(:,i))/l2norm(cross(tV(:,i),cV(:,i)));
end

if circBool
    pV(:,1) = cross(tV(:,1),cV(:,1))/l2norm(cross(tV(:,1),cV(:,1)));
end
%% q vectors
qV = zeros(3,P);
for i = 1:P
    qV(:,i) = cross(tV(:,i),pV(:,i));
end

%% torsion and bond angles

if circBool
    cosPsi = zeros(1,P);
    cosTheta = zeros(1,P);
    cosPsi(P) = dot(tV(:,1), tV(:,P));
    cosTheta(P) = dot(bV(:,1), bV(:,P));
    cosPsi(1) = dot(tV(:,2), tV(:,1));
    cosTheta(1) = dot(bV(:,2), bV(:,1));
    for i = 2:P-1
        cosPsi(i) = dot(tV(:,i+1), tV(:,i));
        cosTheta(i) = dot(bV(:,i+1), bV(:,i));
    end

else
    cosPsi = zeros(1,P-3);
    cosTheta = zeros(1,P-3);
    for i = 2:P-2
        cosPsi(i-1) = dot(tV(:,i+1), tV(:,i));
        cosTheta(i-1) = dot(bV(:,i+1), bV(:,i));
    end
end
if circBool
    cosPsi(P) = dot(tV(:,1), tV(:,P));
    cosTheta(P) = dot(bV(:,1), bV(:,P));
end

bondAngle = acosd(cosPsi);
torsionAngle = acosd(cosTheta);



    function [d] = euclidDist(r2,r1)
        d = sqrt(  (r2(1) - r1(1))^2 + (r2(2) - r1(2))^2 + (r2(3) - r1(3))^2    );
    end

    function [l2] = l2norm(x)
        l2 = abs(sqrt(sum(x.^2)));

    end

%% calculate chi, psi
costc = max(min(dot(tV,cV)/(norm(tV)*norm(cV)),1),-1);
chiAngles = real(acosd(costc));

norm_nb = cross(nV,bV);
proj_cnb = cV - (dot(cV,norm_nb)/norm(norm_nb)).*norm_nb;

cosnc = max(min(dot(nV,proj_cnb)/(norm(nV)*norm(proj_cnb)),1),-1);
psiAngles = real(acosd(cosnc));


%% plot both bases
if allPlotBool
    aS = 0.3; %arrow length
    figure()
    if allPlotBool
        plot3(rAll(:,1), rAll(:,2), rAll(:,3), 'k.','MarkerSize',9)
        hold on
    end
    hold on
    if circBool
        plot3([r(end,1) r(1,1)], [r(end,2) r(1,2)], [r(end,3) r(1,3)], 'Linewidth',5,'Color','k')
    end
    plot3(r(:,1), r(:,2), r(:,3), 'Linewidth',5,'Color','k')

    vwidth = 1.5;
    quiver3(r(:,1),r(:,2),r(:,3),tV(1,:)',tV(2,:)',tV(3,:)',aS,'r','Linewidth',vwidth)
    quiver3(r(:,1),r(:,2),r(:,3),bV(1,:)',bV(2,:)',bV(3,:)',aS,'b','Linewidth',vwidth)
    quiver3(r(:,1),r(:,2),r(:,3),nV(1,:)',nV(2,:)',nV(3,:)',aS,'g','Linewidth',vwidth)
    % quiver3(r(:,1),r(:,2),r(:,3),cV(1,:)',cV(2,:)',cV(3,:)',aS,'r','Linewidth',vwidth)
    % quiver3(r(:,1),r(:,2),r(:,3),pV(1,:)',pV(2,:)',pV(3,:)',aS,'magenta','Linewidth',vwidth)
    % quiver3(r(:,1),r(:,2),r(:,3),qV(1,:)',qV(2,:)',qV(3,:)',aS,'cyan','Linewidth',vwidth)

    view(11,49)
    set(gca,'Color','None')
    axis off
end

%% plot curvature, torsion
figure()
plot(1:length(bondAngle),bondAngle,'b','Linewidth',3, 'Color','#e7298a')
hold on
plot(1:length(torsionAngle),torsionAngle,'r','Linewidth',3,'Color','#2b8cbe')
xlabel('Atom Step')
ylabel('Angle')
set(gca,'Color','None')
legend({'\kappa','\tau'})
ylim([0 180]);
box('on')
if plotBaseNum > 0
    xlim([1, plotBaseNum])
end
set(gca,'FontSize',14)
set(gca,'FontName','Arial')
saveas(gcf,strcat(pwd, '/Figs/', fname, 'CurvTors.fig'));
saveas(gcf,strcat(pwd, '/Figs/', fname, 'CurvTors.png'));

%% curvature/torsion surface
if allPlotBool
    figure()
    Csurf = colorsurf(r,rS,bondAngle,hot,circBool,'Curvature');
    saveas(gcf,strcat(pwd, '/Figs/', fname, 'Cur_surf.fig'));
    saveas(gcf,strcat(pwd, '/Figs/', fname, 'Cur_surf.png'));

    Tsurf = colorsurf(r,rS,torsionAngle,parula,circBool,'Torsion');
    saveas(gcf,strcat(pwd, '/Figs/', fname, 'Tor_surf.fig'));
    saveas(gcf,strcat(pwd, '/Figs/', fname, 'Tor_surf.png'));

end

%% twist writhe
lk_max = 0.5; lk_min = -0.2;
wr_min = -0.1; wr_max = 0.1;
DDNA = 1; %can fix if you want LDNA to have negative or positive link--default is to allow negative link in LDNA

[Wr_global,Lk_global, Wr_x, Lk_x, Lk_all, Wr_all] = Linking2D(r, rS,circBool);
Wr_global = Wr_global.*DDNA; Lk_global = Lk_global.*DDNA;
Wr_x = Wr_x.*DDNA; Lk_x = Lk_x.*DDNA; Wr_all = Wr_all.*DDNA; Lk_all = Lk_all.*DDNA;
Tw_x = Lk_x-Wr_x;
%save structures
linkdata = {Wr_global,Lk_global, Wr_x, Lk_x, Lk_all, Wr_all};
sname = strcat(fname, '_LinkData.mat');
savefile = (sname);
cd LinkData/
save(savefile, 'linkdata','-mat');
cd ..

disp({fname,Lk_global})
disp({fname,Wr_global})

wrsurf = Wrcolorsurf(r,rS,Wr_x,parula,circBool,'Local Writhe',0);
lksurf = Wrcolorsurf(r,rS,Lk_x,jet,circBool,'Local Link',1);


saveas(wrsurf,strcat(pwd, '/Figs/', fname, 'WrSurf.fig'));
saveas(wrsurf,strcat(pwd, '/Figs/', fname, 'WrSurf.png'));
saveas(lksurf,strcat(pwd, '/Figs/', fname, 'LkSurf.fig'));
saveas(lksurf,strcat(pwd, '/Figs/', fname, 'LkSurf.png'));

figure()
hold on
if plotBaseNum > 0
    xlim([1, plotBaseNum])
end
p1 = plot(0*Lk_x(1:end,1),'Color','k');
pa = plot(Wr_x(1:end,1),'Color','#8c6bb1','Linewidth',3,'DisplayName','Wr_{local}');
pb = plot(Tw_x(1:end,1),'Color','#41ae76','Linewidth',3,'DisplayName','Tw_{local}');
pc = plot(Lk_x(1:end,1),'Color','#91003f','Linewidth',3,'DisplayName','Lk_{local}');
xlabel('Frame Step')
ylabel('')
set(gca,'Color','None')
set(gca,'FontSize',14)
set(gca,'FontName','Arial')
hold off
legend([pa(1),pb(1),pc(1)])
ylim([lk_min lk_max]);
box('on')
saveas(gcf,strcat(pwd, '/Figs/', fname, 'LkWrTw.fig'));
saveas(gcf,strcat(pwd, '/Figs/', fname, 'LkWrTw.png'));

%% 2D contour plot of link
figure
colormap('turbo');
contourf(Lk_all);
title('Linking Map');
colorbar(gca,'LineWidth',1)
xlabel('Frame Step')
ylabel('Frame Step')
set(gca,'Color','None')
xlim([1 size(Lk_all,1)]);
ylim([1 size(Lk_all,2)]);
box('on','LineWidth',1,'Color','black')
set(gca,'XGrid','on')
set(gca,'YGrid','on')
set(gca,'GridColor','white')
set(gca,'GridAlpha',0.3)
set(gca,'FontSize',14)
set(gca,'FontName','Arial')
saveas(gcf,strcat(pwd, '/Figs/', fname, '_Lk2D.fig'));
saveas(gcf,strcat(pwd, '/Figs/', fname, '_Lk2D.png'));

%% 2D countour plot of writhe
figure
colormap('bone')
contourf(Wr_all);
title('Writhe Map');
colorbar
xlabel('Frame Step')
ylabel('Frame Step')
set(gca,'Color','None')
xlim([1 size(Wr_all,1)]);
ylim([1 size(Wr_all,2)]);
box('on')
set(gca,'FontSize',14)
set(gca,'FontName','Arial')
saveas(gcf,strcat(pwd, '/Figs/', fname, '_Wr2D.fig'));
saveas(gcf,strcat(pwd, '/Figs/', fname, '_Wr2D.png'));

%% Functions
function [s] = colorsurf(r,rS,c,cmap,circBool,plotName)
    %set colors
    figure()
    set(gcf,'Units',"inches");
    set(gcf,'Position',[3,3,4*560/420,4]);
    if allPlotBool
        plot3(rAll(:,1), rAll(:,2), rAll(:,3), 'k.','MarkerSize',5)
        hold on
    end
    lw = 2;
    maxc = 180;
    cn = (c-min(c))/(maxc - min(c));
    cn = ceil(cn*floor(size(cmap,1)*1));
    cn = max(cn,1);
    if circBool
        cm = repmat(220,[size([rS(:,3);rS(1,3)],1),2,3]);
    else
        cm = repmat(220,[size([rS(:,3);rS(1,3)],1)-1,2,3]);
    end

    for p = 1:length(cn)
        cm(p,1,:) = cmap(cn(p),:);
        cm(p,2,:) = cmap(cn(p),:);
    end
    colormap(cmap)
    cb = colorbar(...'Ticks',0:0.16:1,... 'TickLabels', 0:30:180,...
        'LineWidth',lw,...
        'Position',[0.849,0.185,0.0276,0.6714],...
        'FontWeight','bold','FontSize',12);
    hold on
    set(cb,'Ticks',220*(1:30:180)./180)
    set(cb,'TickLabels',0:30:180)


    if circBool
        s = surf([[r(:,1);r(1,1)], [rS(:,1);rS(1,1)]],...
            [[r(:,2);r(1,2)], [rS(:,2);rS(1,2)]],...
            [[r(:,3);r(1,3)], [rS(:,3);rS(1,3)]],cm,...
            'Linewidth',lw,'FaceAlpha',0.8);
    else
        s = surf([r(:,1), rS(:,1)],[r(:,2), rS(:,2)],[r(:,3), rS(:,3)],...
            cm, 'Linewidth',lw,'FaceAlpha',0.8);
    end

    view(10,69)
    set(gca,'Color','None')
    title(plotName)
    axis off
end

function [s] = Wrcolorsurf(r,rS,c,cmap,circBool,plotName,lk_bool)
    %lkbool is whether this is link (1) or writhe(0)
    %set colors--fix if overflow
lk_max = 0.5; lk_min = -0.25;
wr_min = -0.1; wr_max = 0.1;


    figure()
    set(gcf,'Units',"inches");
    set(gcf,'Position',[3,3,4*560/420,4]);
    if allPlotBool
        plot3(rAll(:,1), rAll(:,2), rAll(:,3), 'k.','MarkerSize',8)
        hold on
    end
    lw = 2;
    if lk_bool
        maxc = lk_max; minc = lk_min;
    else
        maxc = wr_max; minc = wr_min;
    end
    cn = (c-minc)/(maxc - minc);
    cn = ceil(cn*floor(size(cmap,1)));
    cn = max(cn,1);
    if circBool
        cm = repmat(220,[size([rS(:,3);rS(1,3)],1),2,3]);
    else
        cm = repmat(220,[size([rS(:,3);rS(1,3)],1)-1,2,3]);
    end

    for p = 1:length(cn)
        if cn(p) > 256
            cn(p) = 256;
        end
        cm(p,1,:) = cmap(cn(p),:);
        cm(p,2,:) = cmap(cn(p),:);
    end
    colormap(cmap)
    cb = colorbar(...
        'LineWidth',lw,...
        'Position',[0.909,0.182,0.0276,0.671],...
        'FontWeight','bold','FontSize',12);
    hold on
    set(cb,'Ticks',0:22:220)
    set(cb,'TickLabels',[])


    if circBool
        s = surf([[r(:,1);r(1,1)], [rS(:,1);rS(1,1)]],...
            [[r(:,2);r(1,2)], [rS(:,2);rS(1,2)]],...
            [[r(:,3);r(1,3)], [rS(:,3);rS(1,3)]],cm,...
            'Linewidth',lw,'FaceAlpha',0.8);
    else
        s = surf([r(:,1), rS(:,1)],[r(:,2), rS(:,2)],[r(:,3), rS(:,3)],...
            cm, 'Linewidth',lw,'FaceAlpha',0.8);
    end

    view(0,-72)
    set(gca,'Color','None')
    title(plotName)
    axis off
end

function [Wr_global, Lk_global, Wr_x, Lk_x, Lk_all, Wr_all] ...
            = Linking2D(K,L,circBool)
    Wr_global = 0;
    Lk_global = 0;
    Wr_all = zeros(size(K,1),size(L,1));
    Lk_all = zeros(size(K,1),size(L,1));
    Lk_x = zeros(size(L,1),size(L,1));
    Wr_x = zeros(size(L,1),size(L,1));
    if circBool
        endnum = 0;
    else
        endnum = 1;
    end

    for x = 1:size(K,1)-endnum %add last point from xend to 1; x is j
        if x == size(K,1) %%circBool implied by endnum--only true for cbool
            K1 = K(end,:);
            K2 = K(1,:);
        else
            K1 = K(x,:);
            K2 = K(x+1,:);
        end
        for y = 1:size(L,1)-endnum %add last point for yend to 1; y is k
            if y == size(L,1)
                L1 = L(end,:);
                L2 = L(1,:);
            else
                L1 = L(y,:);
                L2 = L(y+1,:);
            end
            dq = L2-L1;
            dp = K2-K1;
            dpxdq = cross(dp,dq);
            rpq = @(t,s) (s*L2+(1-s)*L1)-(t*K2+(1-t)*K1);
            funrpq = @(t,s) dot(dpxdq,rpq(t,s)/(norm(rpq(t,s),2)^3));
            intr = integral(@(t) integral(@(s) funrpq(t,s),0,1,...
                'ArrayValued',true),0,1,'ArrayValued',true)/(4*pi);

            if x ~= y && x ~= y-1 && x ~= y+1
                Wr_global = Wr_global - intr; %% minus sign flipped here
                Wr_all(x,y) = -1*intr;
                Wr_x(x) = Wr_x(x,1) - intr;
            end
            Lk_global = Lk_global - intr; %minus sign flipped here
            Lk_all(x,y) = -1*intr;
            Lk_x(x) = Lk_x(x,1) - intr;
        end

    end
end

end

