function [Patch,Area1,Area2] = PatchGenerate3(ActiveVoxSeed,Patch,VertConn,VertArea,overlap_area,area2)
overlap_area=overlap_area*1e-4;
area2=area2*1e-4;
Area1 = sum(VertArea(Patch));
Area2=0;
while Area1 <= overlap_area && Area1 <= area2
    if Area1 <= overlap_area && Area1 <= area2
        verts = find(max(VertConn(Patch, :), [], 1));
        newverts = setdiff(verts, Patch);%新增点
        over_verts=intersect(ActiveVoxSeed,newverts);%新增重叠点
        if ~isempty(over_verts)
            over_Nouter = union(Patch,over_verts);%重叠点
        else
            break;
        end
        Area1 = sum(VertArea(over_Nouter));
    end
    if Area1 > overlap_area || Area1 > area2
        overlap_Ndiff =  setdiff(over_Nouter,Patch);%新增重叠点
        for i = 1: numel(overlap_Ndiff)
            Patch = union(Patch,overlap_Ndiff(i));
            Area1 = sum(VertArea(Patch));
            if Area1 > overlap_area || Area1 > area2
                break;
            end
        end
    else
        Patch=over_Nouter;
    end
end
Patch1=Patch;
diff_area=area2-Area1;
if diff_area<=0
    return 
end
while Area2 <= diff_area
    if Area2 <= diff_area
        verts = find(max(VertConn(Patch, :), [], 1));
        newverts = setdiff(verts, Patch);%新增点
        diff_verts=setdiff(newverts,ActiveVoxSeed);%新增非重叠点
        if ~isempty(diff_verts)
            Nouter = union(Patch,diff_verts);%重叠点+非重叠点
        else
            break;
        end
        Area2 = sum(VertArea(setdiff(Nouter,Patch1)));
    end
    if Area2 > diff_area
        diff_Ndiff =  diff_verts;%新增非重叠点
        for i = 1: numel(diff_Ndiff)
            Patch = union(Patch,diff_Ndiff(i));
            Area2 =Area2+ sum(VertArea(diff_Ndiff(i)));
            if Area2 > diff_area
                break;
            end
        end
    else
        Patch=Nouter;
    end
end

