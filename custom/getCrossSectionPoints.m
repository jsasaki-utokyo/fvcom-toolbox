function inputConf = getCrossSectionPoints(inputConf)

    if length(inputConf.endpoints)~=length(inputConf.crosssection)
        error('Input cross section number dismatch endpoints, check again!'); 
    end

    ds = inputConf.crosssection_ds;

    course = [];
    name = {};
    num = [];
    for i = 1:length(inputConf.crosssection)
        pointname = inputConf.crosssection{i};
        endpoint1 = [inputConf.endpoints(i,1),inputConf.endpoints(i,2)];
        endpoint2 = [inputConf.endpoints(i,3),inputConf.endpoints(i,4)];
        segments = round(sqrt((endpoint1(1)-endpoint2(1))^2+(endpoint1(2)-endpoint2(2))^2)/ds);
        course_x = linspace(endpoint1(1),endpoint2(1),segments);
        course_y = linspace(endpoint1(2),endpoint2(2),segments);
        for j = 1:segments
            temp1(j,1) = course_x(j);
            temp1(j,2) = course_y(j);
            temp2{j} = [pointname,'_',num2str(j,'%04d')];
        end
        course = [course;temp1];
        name = [name,temp2];
        num = [num,segments];
        clear temp*;
    end

    inputConf.positions = [inputConf.positions; course];
    for k = 1:length(name)
        inputConf.names{end+1} = name{k};
    end
end