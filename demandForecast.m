function result = demandForecast(dbName,dbID,dbPW,VER_ID,...
    periodStart,periodEnd,fcstStart,fcstEnd,varargin)

parseObj = inputParser;
parseObj.addRequired('dbName');
parseObj.addRequired('dbID');
parseObj.addRequired('dbPW');
parseObj.addRequired('VER_ID');
parseObj.addRequired('periodStart');
parseObj.addRequired('periodEnd');
parseObj.addRequired('fcstStart');
parseObj.addRequired('fcstEnd');
parseObj.addParameter('DSTRB_CHNL',"ALL");
parseObj.addParameter('SALES_GRP',"ALL");
parseObj.addParameter('parameterInit',0);
parseObj.addParameter('maxLag',3);
parseObj.addParameter('maxLagS',2);
parseObj.addParameter('numTest',12);
% parseObj.addParameter('numForecasts',3);

parseObj.parse(dbName,dbID,dbPW,VER_ID,...
    periodStart,periodEnd,fcstStart,fcstEnd,varargin{:});

dbName = parseObj.Results.dbName;
dbID = parseObj.Results.dbID;
dbPW = parseObj.Results.dbPW;
VER_ID = parseObj.Results.VER_ID;
periodStart = parseObj.Results.periodStart;
periodEnd = parseObj.Results.periodEnd;
fcstStart = parseObj.Results.fcstStart;
fcstEnd = parseObj.Results.fcstEnd;
DSTRB_CHNL = parseObj.Results.DSTRB_CHNL;
SALES_GRP = parseObj.Results.SALES_GRP;
parameterInit = parseObj.Results.parameterInit;
maxLag = parseObj.Results.maxLag;
maxLagS = parseObj.Results.maxLagS;
numTest = parseObj.Results.numTest;
% numForecasts = parseObj.Results.numForecasts;

[OK01,result01] = dbNameCheck(dbName);
if OK01 == false
    result = result01;
end

[OK02,result02] = dbIDCheck(dbID);
if OK02 == false
    result = result02;
end

[OK03,result03] = dbPWCheck(dbPW);
if OK03 == false
    result = result03;
end

[OK04,result04] = VER_IDCheck(VER_ID);
if OK04 == false
    result = result04;
end

[OK05,result05] = periodCheck(periodStart);
if OK05 == false
    result = result05;
end

[OK06,result06] = periodCheck(periodEnd);
if OK06 == false
    result = result06;
end

[OK07,result07] = periodCheck(fcstStart);
if OK07 == false
    result = result07;
end

[OK08,result08] = periodCheck(fcstEnd);
if OK08 == false
    result = result08;
end

[OK09,result09] = DSTRB_CHNLCheck(DSTRB_CHNL);
if OK09 == false
    result = result09;
end

[OK10,result10] = SALES_GRPCheck(SALES_GRP);
if OK10 == false
    result = result10;
end

[OK11,result11] = parameterInitCheck(parameterInit);
if OK11 == false
    result = result11;
end

[OK12,result12] = maxLagCheck(maxLag);
if OK12 == false
    result = result12;
end

[OK13,result13] = maxLagCheck(maxLagS);
if OK13 == false
    result = result13;
end

[OK14,result14] = numTestCheck(numTest);
if OK14 == false
    result = result14;
end

% [OK15,result15] = numTestCheck(numForecasts);
% if OK15 == false
%     result = result15;
% end


warning('off','all')


if (OK01 * OK02 * OK03 * OK04 * OK05 * OK06 * OK07 * OK08 * OK09 *...
        OK10 * OK11 * OK12 * OK13 * OK14) == 1
    
    conn = database(dbName, dbID, dbPW);

    queryDataBase = strcat("SELECT * FROM fcst.dbo.FCST_ENGINE_ORD_RSLT_INF WHERE VER_ID = '",...
        VER_ID, "'");
    queryParameter = "SELECT * FROM fcst.dbo.FCST_ENGINE_SAVED_PARAMETER";
%     queryForecast = ['SELECT * FROM fcst.dbo.FCST_ENGINE_RSLT'];
    queryMasterBase = strcat("SELECT YYYY FROM fcst.dbo.FCST_VER_MST_DTL WHERE USE_YN = 'N' AND VER_ID = '",...
        VER_ID, "'");

%     tblRaw = fetch(conn,queryData);
    tblParameter = fetch(conn,queryParameter);
%     tblForecast = fetch(conn,queryForecast);

%     targetVer = strcmp(tblRaw.VER_ID,VER_ID);
%     targetVerIndex = find(targetVer);
%     tblVer = tblRaw(targetVerIndex,:);

    if (DSTRB_CHNL == "ALL") && (SALES_GRP == "ALL")
%         tbl = tblVer;
        queryData = queryDataBase;
        queryMaster = queryMasterBase;

    elseif (DSTRB_CHNL == "ALL") && (SALES_GRP ~= "ALL")
%         targetPos = strcmp(tblVer.FCST_SALES_GRP_ID,SALES_GRP);
%         targetIndex = find(targetPos);
%         tbl = tblVer(targetIndex,:);
        queryPlus = strcat(" AND FCST_SALES_GRP_ID = '", SALES_GRP, "'");
        queryData = strcat(queryDataBase, queryPlus);
        queryMaster = strcat(queryMasterBase,queryPlus);

    elseif (DSTRB_CHNL ~= "ALL") && (SALES_GRP == "ALL")
%         targetPos = strcmp(tblVer.FCST_DSTRB_CHNL_ID,DSTRB_CHNL);
%         targetIndex = find(targetPos);
%         tbl = tblVer(targetIndex,:);
        queryPlus = strcat(" AND FCST_DSTRB_CHNL_ID = '", DSTRB_CHNL, "'");
        queryData = strcat(queryDataBase, queryPlus);
        queryMaster = strcat(queryMasterBase,queryPlus);

    else
%         targetPos = strcmp(tblVer.FCST_DSTRB_CHNL_ID,DSTRB_CHNL) .*...
%             strcmp(tblVer.FCST_SALES_GRP_ID,SALES_GRP);
%         targetIndex = find(targetPos);
%         tbl = tblVer(targetIndex,:);
        queryPlus = strcat(" AND FCST_DSTRB_CHNL_ID = '", DSTRB_CHNL, "'",...
            " AND FCST_SALES_GRP_ID = '", SALES_GRP, "'");
        queryData = strcat(queryDataBase, queryPlus);
        queryMaster = strcat(queryMasterBase,queryPlus);

    end
    
    tbl = fetch(conn, queryData);
    tblMaster = fetch(conn, queryMaster);

    grp = findgroups(tbl.FCST_DSTRB_CHNL_ID,...
        tbl.FCST_SALES_GRP_ID, tbl.HIER_LVL_4_ID);
    
    delbase = "DELETE FROM fcst.dbo.FCST_ENGINE_RSLT WHERE VER_ID = '";
    deldstrb = "'and FCST_DSTRB_CHNL_ID = '";
    delsales = "'and FCST_SALES_GRP_ID = '";

    if ~isempty(grp)
        tblSplit = splitapply(@(varargin) varargin, tbl, grp);

        numID = size(tblSplit,1);
        tblSub = cell(numID,1);

        for i = 1:numID
            tblGroup = table(tblSplit{i,:},...
                'VariableNames',tbl.Properties.VariableNames);
            tblGroupSorted = sortrows(tblGroup,'YYYYMM');

            tblSub{i} = tblGroupSorted;
        end

        period = 12;

        % tblParameterNew = [];
        % tblForecastNew = [];
        
        if (DSTRB_CHNL == "ALL") && (SALES_GRP == "ALL")
            delquery = strcat(delbase, VER_ID, "'");

        elseif (DSTRB_CHNL == "ALL") && (SALES_GRP ~= "ALL")
            delquery = strcat(delbase, VER_ID, delsales, SALES_GRP, "'");

        elseif (DSTRB_CHNL ~= "ALL") && (SALES_GRP == "ALL")
            delquery = strcat(delbase, VER_ID, deldstrb, DSTRB_CHNL, "'");

        else
            delquery = strcat(delbase, VER_ID, deldstrb, DSTRB_CHNL, delsales, SALES_GRP, "'");
            
        end
        
        execute(conn, delquery);

        for i = 1:numID
            
            periodStartDT = datetime(periodStart,'InputFormat','yyyyMM');
            periodEndDT = datetime(periodEnd,'InputFormat','yyyyMM');
            fcstStartDT = datetime(fcstStart,'InputFormat','yyyyMM');
            fcstEndDT = datetime(fcstEnd,'InputFormat','yyyyMM');
            
            periodLen = calmonths(between(periodStartDT,periodEndDT));
            btwForecasts = calmonths(between(fcstStartDT,fcstEndDT));
            
            numForecasts = calmonths(between(periodEndDT,fcstEndDT));
            fcstFirstNum = numForecasts - btwForecasts;
            
            numForecastsReport = btwForecasts + 1;
            
            wholeMonths = periodStartDT + calmonths(0:periodLen);
            fcstMonths = fcstStartDT + calmonths(0:btwForecasts);

            sample = tblSub{i};
            sampleNeed = sample(:,{'YYYYMM','ORD_QTY'});
            
            sampleQTY = zeros((periodLen+1),1);
            sampleDelIdx = zeros((periodLen+1),1);
            
            yearDel = tblMaster.YYYY;
            
            for j = 1:(periodLen+1)
                periodDate = wholeMonths(j);
                sampleYear = num2str(year(periodDate));
                
                if any(strcmp(yearDel,sampleYear))
                    sampleDelIdx(j,1) = 1;
                else
                    sampleMonth = num2str(month(periodDate),'%02.f');
                    sampleDate = strcat(sampleYear, sampleMonth);

                    monthPos = strcmp(sampleNeed.YYYYMM,string(sampleDate));
                    monthIdx = find(monthPos);

                    if ~isempty(monthIdx)
                        sampleQTY(j,1) = sampleNeed.ORD_QTY(monthIdx(1));
                    end
                end
            end
            
            sampleQTY(find(sampleDelIdx),:) = [];
            data = fillmissing(sampleQTY,'constant',0);
            
%             verID = VER_ID;
            verID = strcat(num2str(year(datetime)),...
                num2str(month(datetime),'%02.f'));
            chnID = char(table2cell(sample(1,'FCST_DSTRB_CHNL_ID')));
            grpID = char(table2cell(sample(1,'FCST_SALES_GRP_ID')));
            tar = char(table2cell(sample(1,'HIER_LVL_4_ID')));
            
            numFill = maxLagS * period + numTest - size(data,1);
            
            if numFill > 0
                for j = 1:numFill
                    
                    targetDate = periodStartDT - calmonths(j);
                    findYear = num2str(year(targetDate));
                    findMonth = num2str(month(targetDate),'%02.f');
                    findDate = strcat(findYear, findMonth);
                    
                    findPos = strcmp(sampleNeed.YYYYMM,findDate);
                    findIdx = find(findPos);

                    if ~isempty(findIdx)
                        data = [sampleNeed.ORD_QTY(findIdx(1));data];
                    
                    else
                        data = [0;data];
                    end
                end
            end
            
            dataInput = data;
            
%             nz = find(data);
%             if isempty(nz)
%                 dataInput = zeros((2*period + numTest),1);
%             elseif (size(data,1)-nz(1)+1) < 2*period + numTest
%                 numzeros = 2*period + numTest - (size(data,1)-nz(1)+1);
%                 dataInput = [zeros(numzeros,1);data(nz(1):end)];
%             else
%                 dataInput = data(nz(1):end);
%             end

            try
                if parameterInit
                    [modelName,parameters,rmse] = modelSelection(dataInput,...
                        period,maxLag,maxLagS,numTest);

                    modelName = convertStringsToChars(modelName);
                    paramVal = num2str(parameters(1));
                    for j = 2:6
                        paramVal = strcat(paramVal, num2str(parameters(j)));
                    end

                    tblParameterTemp(1,'VER_ID') = table({verID});
                    tblParameterTemp(1,'FCST_DSTRB_CHNL_ID') = table({chnID});
                    tblParameterTemp(1,'FCST_SALES_GRP_ID') = table({grpID});
                    tblParameterTemp(1,'HIER_LVL_4_ID') = table({tar});
                    tblParameterTemp(1,'FCST_MODEL') = table({modelName});
                    tblParameterTemp(1,'FCST_PARAMETER') = table({paramVal});
                    tblParameterTemp(1,'FCST_RMSE') = table(rmse);
                    currentTime = strcat(char(datetime),'.0000000');
                    tblParameterTemp(1,'INS_DT') = table({currentTime});
                    tblParameterTemp(1,'INSP_ID') = table({blanks(0)});
                    tblParameterTemp(1,'UPD_DT') = table({currentTime});
                    tblParameterTemp(1,'UPDP_ID') = table({blanks(0)});

        %             tblParameterNew = [tblParameterNew;tblParameterTemp];

                    if ~any(strcmp(tblParameter.FCST_DSTRB_CHNL_ID,chnID) .*...
                            strcmp(tblParameter.FCST_SALES_GRP_ID,grpID) .*...
                            strcmp(tblParameter.HIER_LVL_4_ID,tar))

                        sqlwrite(conn,'FCST_ENGINE_SAVED_PARAMETER',tblParameterTemp,...
                            'Catalog','fcst','Schema','dbo');

                    else
                        colParam = {'VER_ID','FCST_MODEL','FCST_PARAMETER',...
                            'FCST_RMSE','UPD_DT','UPDP_ID'};
                        whereParam = strcat("WHERE FCST_DSTRB_CHNL_ID = '",chnID,...
                            "' AND FCST_SALES_GRP_ID = '",grpID,...
                            "' AND HIER_LVL_4_ID = '",tar,...
                            "'");
                        updateParam = tblParameterTemp(1,colParam);

                        update(conn,'fcst.dbo.FCST_ENGINE_SAVED_PARAMETER',colParam,...
                            updateParam,whereParam);
                    end

                else
                    if ~any(strcmp(tblParameter.FCST_DSTRB_CHNL_ID,chnID) .*...
                            strcmp(tblParameter.FCST_SALES_GRP_ID,grpID) .*...
                            strcmp(tblParameter.HIER_LVL_4_ID,tar))

                        [modelName,parameters,rmse] = modelSelection(dataInput,...
                            period,maxLag,maxLagS,numTest);

                        modelName = convertStringsToChars(modelName);
                        paramVal = num2str(parameters(1));
                        for j = 2:6
                            paramVal = strcat(paramVal, num2str(parameters(j)));
                        end

                        tblParameterTemp(1,'VER_ID') = table({verID});
                        tblParameterTemp(1,'FCST_DSTRB_CHNL_ID') = table({chnID});
                        tblParameterTemp(1,'FCST_SALES_GRP_ID') = table({grpID});
                        tblParameterTemp(1,'HIER_LVL_4_ID') = table({tar});
                        tblParameterTemp(1,'FCST_MODEL') = table({modelName});
                        tblParameterTemp(1,'FCST_PARAMETER') = table({paramVal});
                        tblParameterTemp(1,'FCST_RMSE') = table(rmse);
                        currentTime = strcat(char(datetime),'.0000000');
                        tblParameterTemp(1,'INS_DT') = table({currentTime});
                        tblParameterTemp(1,'INSP_ID') = table({blanks(0)});
                        tblParameterTemp(1,'UPD_DT') = table({currentTime});
                        tblParameterTemp(1,'UPDP_ID') = table({blanks(0)});

        %                 tblParameterNew = [tblParameterNew;tblParameterTemp];

                        sqlwrite(conn,'FCST_ENGINE_SAVED_PARAMETER',tblParameterTemp,...
                            'Catalog','fcst','Schema','dbo');

                    else
                        parameterPos = strcmp(tblParameter.FCST_DSTRB_CHNL_ID,chnID) .*...
                            strcmp(tblParameter.FCST_SALES_GRP_ID,grpID) .*...
                            strcmp(tblParameter.HIER_LVL_4_ID,tar);
                        parameterIndex = find(parameterPos);

                        tblParameterSample = tblParameter(parameterIndex,:);
                        tblParameterSampleSorted = sortrows(tblParameterSample,'VER_ID');

                        modelName = tblParameterSampleSorted.FCST_MODEL(end);
                        PARAM = char(tblParameterSampleSorted.FCST_PARAMETER(end));

                        P = str2double(PARAM(1));
                        D = str2double(PARAM(2));
                        Q = str2double(PARAM(3));
                        PS = str2double(PARAM(4));
                        DS = str2double(PARAM(5));
                        QS = str2double(PARAM(6));
                        parameters = [P, D, Q, PS, DS, QS];
                    end
                end

                [~,outForecasts,~,~,~] = forecastFuture(dataInput,...
                    modelName,parameters,period,numForecasts);

            catch
                outForecasts = ones(numForecasts,1) .* (-99999999);
            end

            forecastVer = cell(numForecastsReport,1);
            forecastChn = cell(numForecastsReport,1);
            forecastGrp = cell(numForecastsReport,1);
            forecastID = cell(numForecastsReport,1);
            forecastDates = cell(numForecastsReport,1);
            forecastIDT = cell(numForecastsReport,1);
            forecastIID = cell(numForecastsReport,1);
            forecastUDT = cell(numForecastsReport,1);
            forecastUID = cell(numForecastsReport,1);
            
            for j = 1:numForecastsReport
%                 forecastVer{j} = verID;
                forecastVer{j} = VER_ID;
                forecastChn{j} = chnID;
                forecastGrp{j} = grpID;
                forecastID{j} = tar;

                yearTemp = num2str(year(fcstMonths(j)));
                monthTemp = num2str(month(fcstMonths(j)),'%02.f');
                forecastDates{j} = strcat(yearTemp, monthTemp);
                
                currentTime = strcat(char(datetime),'.0000000');
                forecastIDT{j} = currentTime;
                forecastIID{j} = blanks(0);
                forecastUDT{j} = currentTime;
                forecastUID{j} = blanks(0);
            end
            
            outForecastsReport = outForecasts(fcstFirstNum:numForecasts);

            tblForecastTemp(1:numForecastsReport,'VER_ID') = table(forecastVer);
            tblForecastTemp(1:numForecastsReport,'FCST_DSTRB_CHNL_ID') = table(forecastChn);
            tblForecastTemp(1:numForecastsReport,'FCST_SALES_GRP_ID') = table(forecastGrp);
            tblForecastTemp(1:numForecastsReport,'HIER_LVL_4_ID') = table(forecastID);
            tblForecastTemp(1:numForecastsReport,'YYYYMM') = table(forecastDates);
            tblForecastTemp(1:numForecastsReport,'FCST_QTY') = table(outForecastsReport);
            tblForecastTemp(1:numForecastsReport,'INS_DT') = table(cellstr(forecastIDT));
            tblForecastTemp(1:numForecastsReport,'INSP_ID') = table(forecastIID);
            tblForecastTemp(1:numForecastsReport,'UPD_DT') = table(cellstr(forecastUDT));
            tblForecastTemp(1:numForecastsReport,'UPDP_ID') = table(forecastUID);
            
        %     tblForecastNew = [tblForecastNew;tblForecastTemp];

            for j = 1:numForecastsReport

%                 if ~any(strcmp(tblForecast.VER_ID,VER_ID) .*...
%                         strcmp(tblForecast.FCST_DSTRB_CHNL_ID,chnID) .*...
%                         strcmp(tblForecast.FCST_SALES_GRP_ID,grpID) .*...
%                         strcmp(tblForecast.HIER_LVL_4_ID,tar) .*...
%                         strcmp(tblForecast.YYYYMM,char(forecastDates(j))))

                    sqlwrite(conn,'FCST_ENGINE_RSLT',tblForecastTemp(j,:),...
                        'Catalog','fcst','Schema','dbo');

%                 else
%                     colFCST = {'FCST_QTY','UPD_DT','UPDP_ID'};
%                     whereFCST = strcat("WHERE VER_ID = '",VER_ID,...
%                         "' AND FCST_DSTRB_CHNL_ID = '",chnID,...
%                         "' AND FCST_SALES_GRP_ID = '",grpID,...
%                         "' AND HIER_LVL_4_ID = '",tar,...
%                         "' AND YYYYMM = '",forecastDates(j),...
%                         "'");
%                     updateFCST = tblForecastTemp(j,colFCST);
% 
%                     update(conn,'fcst.dbo.FCST_ENGINE_RSLT',colFCST,...
%                         updateFCST,whereFCST);
%                 end
            end

        end

        s = struct;
        s.status = "success";
        s.message = sprintf('%d건', numID);

        result = jsonencode(s);

    %     fid = fopen('result.json','w');
    %     fprintf(fid,'%s',result);
    %     fclose(fid);

    else
        s = struct;
        s.status = "error";
        s.message = "No row for such key is found in the table.";

        result = jsonencode(s);

    %     fid = fopen('result.json','w');
    %     fprintf(fid,'%s',result);
    %     fclose(fid);
    end

    close(conn)
    clear conn query
end


function [adfstat,pval,critval,resid,lags,ICs]=augdfautolag(y,p,maxlags) % ,IC)

% Setup common to all problems
%y=y-y(1);
ydiff=diff(y);
[ydiffcurr, ydifflags]=newlagmatrix(ydiff,maxlags); %#ok<ASGLU>
T=length(y);
Y=y(maxlags+2:T);
tau=length(Y);
%
switch p
    case 0
        % Case 1
        i=0;
        X=y(maxlags+1:T-1);
        rho = X\Y;
        % Compute the errors
        e= Y-X*rho;
        s2(i+1)=e'*e/tau;
        K(i+1)=size(X,2);
        
        % Loop
        for i=1:maxlags
            X=[y(maxlags+1:T-1) ydifflags(:,1:i)];
            rho = X\Y;
            % Compute the errors
            e= Y-X*rho;
            s2(i+1)=e'*e/tau;
            K(i+1)=size(X,2);
        end
        
    case {1,3}
        %Case 2
        i=0;
        X=[ones(size(Y)) y(maxlags+1:T-1)];
        rho = X\Y;
        % Compute the errors
        e= Y-X*rho;
        s2(i+1)=e'*e/tau;
        K(i+1)=size(X,2);
        
        % Loop
        for i=1:maxlags
            X=[ones(size(Y)) y(maxlags+1:T-1) ydifflags(:,1:i)];
            rho = X\Y;
            % Compute the errors
            e= Y-X*rho;
            s2(i+1)=e'*e/tau;
            K(i+1)=size(X,2);
        end
        
        
    case 2
        %Case 4
        i=0;
        X=[ones(size(Y)) y(maxlags+1:T-1) (1:tau)'];
        rho = X\Y;
        % Compute the errors
        e= Y-X*rho;
        s2(i+1)=e'*e/tau;
        K(i+1)=size(X,2);
        
        % Loop
        for i=1:maxlags
            X=[ones(size(Y)) y(maxlags+1:T-1) (1:tau)' ydifflags(:,1:i)];
            rho = X\Y;
            % Compute the errors
            e= Y-X*rho;
            s2(i+1)=e'*e/tau;
            K(i+1)=size(X,2);
        end
end


% if strcmp(IC,'AIC')
    ICs=log(s2) + 2*K/tau;
% else
%     ICs=log(s2) + K*log(tau)/tau;
% end
[~,lags]=min(ICs);
lags=lags-1;
[adfstat,pval,critval,resid]=augdf(y,p,lags);


function [y,x]=newlagmatrix(x,nlags) % ,c)

T=size(x,1);

if nlags>0
    nlags=nlags+1;
    newX=[x;zeros(nlags,1)];
    lagmatrix=repmat(newX,nlags,1);
    lagmatrix=reshape(lagmatrix(1:size(lagmatrix,1)-nlags),T+nlags-1,nlags);
    lagmatrix=lagmatrix(nlags:T,:);
    y=lagmatrix(:,1);
    x=lagmatrix(:,2:nlags);
%     if c==1
%         x=[ones(size(x,1),1) x];
%     end
else
%     if c==1
%         y=x;
%         x=ones(T,1);
%     else
        y=x;
        x=[];
%     end
end


function [adfstat,pval,critval,resid]=augdf(y,p,lags)

% Setup common to all problems
%y=y-y(1);
ydiff=diff(y);
[ydiffcurr, ydifflags]=newlagmatrix(ydiff,lags);
T=length(y);
Y=y(lags+2:T);
tau=length(Y);
%
switch p
    case 0
        % Case 1
        X=[y(lags+1:T-1) ydifflags];
        rho = X\ydiffcurr;
        % Compute the errors
        e= ydiffcurr-X*rho;

        % Compute the covariance matrix of the estimated parameters
        s2 = e'*e/(tau-size(X,2));
        Uinv=inv(diag([T T^(0.5)*ones(1,lags)]));
        sel=[1 zeros(1,lags)];
        sigp=s2*sel*(Uinv*(X'*X)*Uinv)^(-1)*sel';

        % Compute the ADF
        adfstat = tau*rho(1)/sqrt(sigp);
        % Look up the pval and critical values
        [pval,critval]=augdfcv(adfstat,p,tau);
    case 1
        %Case 2
        X=[ones(size(Y)) y(lags+1:T-1) ydifflags];
        rho = X\ydiffcurr;
        % Compute the errors
        e= ydiffcurr-X*rho;

        % Compute the covariance matrix of the estimated parameters
        s2 = e'*e/(tau-size(X,2));
        Uinv=inv(diag([T^(0.5) T T^(0.5)*ones(1,lags)]));
        sel=[0 1 zeros(1,lags)];
        sigp=s2*sel*(Uinv*(X'*X)*Uinv)^(-1)*sel';

        % Compute the ADF
        adfstat = tau*rho(2)/sqrt(sigp);
        % Look up the pval and critical values
        [pval,critval]=augdfcv(adfstat,p,tau);

    case 2
        %Case 4
        X=[ones(size(Y)) y(lags+1:T-1) (1:tau)' ydifflags];
        rho = X\Y;
        % Compute the errors
        e= Y-X*rho;

        % Compute the covariance matrix of the estimated parameters
        s2 = e'*e/(tau-size(X,2));
        X(:,2)=X(:,2)-rho(1)*(1:tau)';
        Uinv=inv(diag([tau^(0.5) tau tau^(1.5) tau^(0.5)*ones(1,lags)]));
        sel=[0 1 0 zeros(1,lags)];
        sigp=s2*sel*(Uinv*(X'*X)*Uinv)^(-1)*sel';

        % Compute the ADF
        adfstat = tau*(rho(2)-1)/sqrt(sigp);
        % Look up the pval and critical values
        [pval,critval]=augdfcv(adfstat,p,tau);

    case 3
        %Case 3
        X=[ones(size(Y)) y(lags+1:T-1) ydifflags];
        rho = X\Y;
        % Compute the errors
        e= Y-X*rho;

        % Compute the covariance matrix of the estimated parameters
        s2 = e'*e/(tau-size(X,2));
        Uinv=inv(diag([tau^(0.5) tau^(1.5)  tau^(0.5)*ones(1,lags)]));
        sel=[0 1 zeros(1,lags)];
        sigp=s2*sel*(Uinv*(X'*X)*Uinv)^(-1)*sel';

        % Compute the ADF
        adfstat = tau^(1.5)*(rho(2)-1)/sqrt(sigp);
        % Look up the pval and critical values
        critval=norminv([.01 .05 .1 .9 .95 .99]');
        pval = normcdf(adfstat);
end
resid=e;


function [modelName,parameters,rmse] = modelSelection(History,period,maxLag,maxLagS,numTest)

len = length(History);

Train = History(1:(len-numTest));
Test = History((len-numTest+1):len);

% Holt-Winter
[alpha,beta,gamma] = parameterHW(Train,period);
[~,forecastsOutHW,~,~,~]...
    = estimateHW(Train,alpha,beta,gamma,period,numTest);

sseHW = 0;
for i = (1:numTest)
    sseHW = sseHW + ((Test(i) - forecastsOutHW(i))^2);
end
rmseHW = sqrt(sseHW/numTest);


% Unit Root Test
% [adfH,~,~,~] = adftest(Train,'alpha',0.05);
[~,pval,~,~,~,~] = augdfautolag(Train,3,3);
adfH = (pval < 0.05);


% ARIMA
try
    [pA,dA,qA] = parameterARIMA(Train,adfH,maxLag);
    [~,forecastsOutA,~,~,~,~]...
        = estimateARIMA(Train,pA,dA,qA,numTest);
    
    sseARIMA = 0;
    for i = (1:numTest)
        sseARIMA = sseARIMA + ((Test(i) - forecastsOutA(i))^2);
    end
    rmseARIMA = sqrt(sseARIMA/numTest);

catch
    pA = 0;
    dA = 0;
    qA = 0;
    rmseARIMA = inf;
end


% SARIMA
try
    [pS,dS,qS,P,D,Q] = parameterSARIMA(Train,adfH,period,maxLag,maxLagS);
    [~,forecastsOutS,~,~,~,~]...
        = estimateSARIMA(Train,pS,dS,qS,P,D,Q,period,numTest);
    
    sseSARIMA = 0;
    for i = (1:numTest)
        sseSARIMA = sseSARIMA + ((Test(i) - forecastsOutS(i))^2);
    end
    rmseSARIMA = sqrt(sseSARIMA/numTest);

catch
    pS = 0;
    dS = 0;
    qS = 0;
    P = 0;
    D = 0;
    Q = 0;
    rmseSARIMA = inf;
end

RMSE = [rmseHW, rmseARIMA, rmseSARIMA];

if min(RMSE) == rmseHW
    modelName = "HW";
    parameters = [0,0,0,0,0,0];
    rmse = rmseHW;
    
elseif min(RMSE) == rmseARIMA
    modelName = "ARIMA";
    parameters = [pA,dA,qA,0,0,0];
    rmse = rmseARIMA;
    
elseif min(RMSE) == rmseSARIMA
    modelName = "SARIMA";
    parameters = [pS,dS,qS,P,D,Q];
    rmse = rmseSARIMA;
    
end


function [alpha,beta,gamma]...
    = parameterHW(History,period)

param0 = [0.2, 0.2, 0.15]; % parameter = [alpha, beta, gamma]
options = optimoptions(@fmincon,'MaxFunEvals',1e100, ...
  'MaxIter',1e100,'TolX',1e-100,'TolFun',1e-100,...
  'OptimalityTolerance',1e-100,'Display','none');
[param,~] = fmincon(@(param) sseHW(History,period,param),...
    param0,[],[],[],[],[0,0,0],[1,1,1],[],options);

alpha = param(1);
beta = param(2);
gamma = param(3);


function SSE = sseHW(data,period,param) % param = [alpha,beta,gamma]
len = length(data);

[forecastsIn,~,~,~,~]...
    = estimateHW(data,param(1),param(2),param(3),period,0); % yhat

SSE = 0;
for i=1:len
    SSE = SSE + ((data(i) - forecastsIn(i)) ^ 2); % {y(i)-yhat(i)}^2
end


function [p,d,q] = parameterARIMA(History,adfH,maxLag)

if adfH == 0
    d = 1;
else
    d = 0;
end

[arlags,malags] = findARIMA(History,d,maxLag);

p = max(arlags);
q = max(malags);


function [arlags,malags] = findARIMA(data,d,maxLag)
len = length(data);
% candidateP = 1:period;
% candidateQ = 1:period;

comb = cell((maxLag+1),1);
comb{1} = 0;

for i = 2:(maxLag+1)
%     candidates = nchoosek(1:maxLag,(i-1));
%     comb{i} = candidates;
    
    % only consider the number of lags
    comb{i} = 1:(i-1);
end

aicScore = inf;
options = optimoptions(@fmincon,'MaxIterations',10);

for i = 1:(maxLag+1)
    numCandidatesP = size(comb{i},1);
    for j = 1:(maxLag+1)
        numCandidatesQ = size(comb{j},1);
        
        for m = 1:numCandidatesP
            testP = comb{i}(m,:);
            
            for n = 1:numCandidatesQ
                testQ = comb{j}(n,:);
                
                if (testP(1) == 0) && (testQ(1) == 0)
                    continue;
                else
                    if (testP(1) == 0)
                        Mdl = arima('D',d,'MALags',testQ);
                    elseif (testQ(1) == 0)
                        Mdl = arima('D',d,'ARLags',testP);
                    else
                        Mdl = arima('D',d,'ARLags',testP,'MALags',testQ);
                    end
                end
                
                % % If heteroscedasticity assumed
                % CondVarMdl = garch(1,1);
                % Mdl.Variance = CondVarMdl;

                try % for several P, Q, error occurs
                    [~,~,logL,~] = estimate(Mdl,data,...
                        'Options',options,'Display','off');
                    
                    numParam = (testP(1)~=0)*size(testP,2) + (testQ(1)~=0)*size(testQ,2) + 1;
                    aic = aicbic(logL,numParam,len);

                    if aic < aicScore
                        aicScore = aic;
                        arlags = testP;
                        malags = testQ;
                    end
                
                catch
                    continue
                end
            end
        end
    end
end


function [p,d,q,P,D,Q]...
    = parameterSARIMA(History,adfH,period,maxLag,maxLagS)

if adfH == 0
    HistorySD = History((period+1):end) - History(1:(end-period));
%     [adfHSD,~,~,~] = adftest(HistorySD,'alpha',0.05);
    [~,pvalSD,~,~,~,~] = augdfautolag(HistorySD,3,3);
    adfHSD = (pvalSD < 0.05);

    if adfHSD == 0
%         [adfHD,~,~,~] = adftest(diff(History),'alpha',0.05);
        [~,pvalD,~,~,~,~] = augdfautolag(diff(History),3,3);
        adfHD = (pvalD < 0.05);

        if adfHD == 0
            d = 1; D = 1;
        else
            d = 1; D = 0;
        end

    else
        d = 0; D = 1;
    end
else
    d = 0; D = 0;
end

[arlags,sarlags,malags,smalags,~] = findSARIMA(History,d,maxLag,maxLagS,period,D);

p = max(arlags);
q = max(malags);
P = max(sarlags);
Q = max(smalags);


function [arlags,sarlags,malags,smalags,aicScore]...
    = findSARIMA(data,d,maxLag,maxLagS,period,seasonDiff)
len = length(data);
% candidateP = 1:period;
% candidateQ = 1:period;

comb = cell((maxLag+1),1);
comb{1} = 0;
for i = 2:(maxLag+1)
%     candidates = nchoosek(1:maxLag,(i-1));
%     comb{i} = candidates;
    
    % only consider the number of lags
    comb{i} = 1:(i-1);
end

combS = cell((maxLagS+1),1);
combS{1} = 0;
for i = 2:(maxLagS+1)
    combS{i} = 1:(i-1);
end

aicScore = inf;
options = optimoptions(@fmincon,'MaxIterations',10);

for i = 1:(maxLag+1)
    numCandidatesP = size(comb{i},1);
    
    for j = 1:(maxLag+1)
        numCandidatesQ = size(comb{j},1);
        
        for m = 1:numCandidatesP
            testP = comb{i}(m,:);
            
            for n = 1:numCandidatesQ
                testQ = comb{j}(n,:);
                
                if (testP(1) == 0) && (testQ(1) == 0)
                    
                    for iS = 1:(maxLagS+1)
                        testPS = combS{iS};

                        for jS = 1:(maxLagS+1)
                            testQS = combS{jS};
                            
                            if (testPS(1) ~= 0) || (testQS(1) ~= 0)

                                if (testPS(1) == 0)
                                    % Following MATLAB notation
                                    Mdl = arima('D',d,'SMALags',period*testQS);

                                elseif (testQS(1) == 0)
                                    % Following MATLAB notation
                                    Mdl = arima('D',d,'SARLags',period*testPS);
                                else
                                    % Following MATLAB notation
                                    Mdl = arima('D',d,'SARLags',period*testPS,'SMALags',period*testQS);
                                end
                                
                                if seasonDiff == 1
                                    Mdl.Seasonality = period;
                                end
                            
                            else
                                if seasonDiff == 1
                                    Mdl = arima('D',d);
                                    Mdl.Seasonality = period;
                                else
                                    continue;
                                end
                            end
                            
                            % % If heteroscedasticity assumed
                            % CondVarMdl = garch(1,1);
                            % Mdl.Variance = CondVarMdl;

                            try % for several P, Q, error occurs
                                [~,~,logL,~] = estimate(Mdl,data,...
                                    'Options',options,'Display','off');
                                
                                numParam = (testP(1)~=0)*size(testP,2) + (testQ(1)~=0)*size(testQ,2)...
                                    + (testPS(1)~=0)*size(testPS,2) + (testQS(1)~=0)*size(testQS,2) + 1;
                                aic = aicbic(logL,numParam,len);

                                if aic < aicScore
                                    aicScore = aic;

                                    arlags = testP;
                                    malags = testQ;

                                    sarlags = testPS;
                                    smalags = testQS;
                                end
                            
                            catch
                                continue;
                            end
                        end
                    end
                else
                    if (testP(1) == 0)
                        for iS = 1:(maxLagS+1)
                            testPS = combS{iS};

                            for jS = 1:(maxLagS+1)
                                testQS = combS{jS};

                                if (testPS(1) ~= 0) || (testQS(1) ~= 0)

                                    if (testPS(1) == 0)
                                        % Following MATLAB notation
                                        Mdl = arima('D',d,'MALags',testQ,'SMALags',period*testQS);

                                    elseif (testQS(1) == 0)
                                        % Following MATLAB notation
                                        Mdl = arima('D',d,'MALags',testQ,'SARLags',period*testPS);
                                    else
                                        % Following MATLAB notation
                                        Mdl = arima('D',d,'MALags',testQ,'SARLags',period*testPS,'SMALags',period*testQS);
                                    end
                                    
                                    if seasonDiff == 1
                                        Mdl.Seasonality = period;
                                    end
                                
                                else
                                    if seasonDiff == 1
                                        Mdl = arima('D',d,'MALag',testQ);
                                        Mdl.Seasonality = period;
                                    else
                                        continue;
                                    end
                                end

                                % % If heteroscedasticity assumed
                                % CondVarMdl = garch(1,1);
                                % Mdl.Variance = CondVarMdl;

                                try % for several P, Q, error occurs
                                    [~,~,logL,~] = estimate(Mdl,data,...
                                        'Options',options,'Display','off');
                                    
                                numParam = (testP(1)~=0)*size(testP,2) + (testQ(1)~=0)*size(testQ,2)...
                                    + (testPS(1)~=0)*size(testPS,2) + (testQS(1)~=0)*size(testQS,2) + 1;
                                aic = aicbic(logL,numParam,len);

                                if aic < aicScore
                                    aicScore = aic;

                                    arlags = testP;
                                    malags = testQ;

                                    sarlags = testPS;
                                    smalags = testQS;
                                end
                                
                                catch
                                    continue;
                                end
                            end
                        end
                    elseif (testQ(1) == 0)
                        for iS = 1:(maxLagS+1)
                            testPS = combS{iS};

                            for jS = 1:(maxLagS+1)
                                testQS = combS{jS};

                                if (testPS(1) ~= 0) || (testQS(1) ~= 0)

                                    if (testPS(1) == 0)
                                        % Following MATLAB notation
                                        Mdl = arima('ARLags',testP,'D',d,'SMALags',period*testQS);

                                    elseif (testQS(1) == 0)
                                        % Following MATLAB notation
                                        Mdl = arima('ARLags',testP,'D',d,'SARLags',period*testPS);
                                    else
                                        % Following MATLAB notation
                                        Mdl = arima('ARLags',testP,'D',d,'SARLags',period*testPS,'SMALags',period*testQS);
                                    end
                                    
                                    if seasonDiff == 1
                                        Mdl.Seasonality = period;
                                    end
                                
                                else
                                    if seasonDiff == 1
                                        Mdl = arima('ARLags',testP,'D',d);
                                        Mdl.Seasonality = period;
                                    else
                                        continue;
                                    end
                                end
                                
                                % % If heteroscedasticity assumed
                                % CondVarMdl = garch(1,1);
                                % Mdl.Variance = CondVarMdl;

                                try % for several P, Q, error occurs
                                    [~,~,logL,~] = estimate(Mdl,data,...
                                        'Options',options,'Display','off');
                                    
                                    numParam = (testP(1)~=0)*size(testP,2) + (testQ(1)~=0)*size(testQ,2)...
                                        + (testPS(1)~=0)*size(testPS,2) + (testQS(1)~=0)*size(testQS,2) + 1;
                                    aic = aicbic(logL,numParam,len);

                                    if aic < aicScore
                                        aicScore = aic;

                                        arlags = testP;
                                        malags = testQ;

                                        sarlags = testPS;
                                        smalags = testQS;
                                    end
                                
                                catch
                                    continue
                                end
                            end
                        end
                    else
                        for iS = 1:(maxLagS+1)
                            testPS = combS{iS};

                            for jS = 1:(maxLagS+1)
                                testQS = combS{jS};

                                if (testPS(1) ~= 0) || (testQS(1) ~= 0)

                                    if (testPS(1) == 0)
                                        % Following MATLAB notation
                                        Mdl = arima('ARLags',testP,'D',d,'MALags',testQ,'SMALags',period*testQS);

                                    elseif (testQS(1) == 0)
                                        % Following MATLAB notation
                                        Mdl = arima('ARLags',testP,'D',d,'MALags',testQ,'SARLags',period*testPS);
                                    else
                                        % Following MATLAB notation
                                        Mdl = arima('ARLags',testP,'D',d,'MALags',testQ,'SARLags',period*testPS,'SMALags',period*testQS);
                                    end
                                    
                                    if seasonDiff == 1
                                        Mdl.Seasonality = period;
                                    end
                                
                                else
                                    if seasonDiff == 1
                                        Mdl = arima('ARLags',testP,'D',d,'MALags',testQ);
                                        Mdl.Seasonality = period;
                                    else
                                        continue;
                                    end
                                end
                                
                                % % If heteroscedasticity assumed
                                % CondVarMdl = garch(1,1);
                                % Mdl.Variance = CondVarMdl;

                                try % for several P, Q, error occurs
                                    [~,~,logL,~] = estimate(Mdl,data,...
                                        'Options',options,'Display','off');
                                    
                                    numParam = (testP(1)~=0)*size(testP,2) + (testQ(1)~=0)*size(testQ,2)...
                                        + (testPS(1)~=0)*size(testPS,2) + (testQS(1)~=0)*size(testQS,2) + 1;
                                    aic = aicbic(logL,numParam,len);

                                    if aic < aicScore
                                        aicScore = aic;

                                        arlags = testP;
                                        malags = testQ;

                                        sarlags = testPS;
                                        smalags = testQS;
                                    end
                                    
                                catch
                                    continue
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


function [forecastsIn,forecastsOut,levelES,trendES,seasonES]...
    = estimateHW(data,alpha,beta,gamma,period,numForecasts)

len = length(data);

% trend0 = (sum(data((period+1):(period+period))) - sum(data(1:period)))...
%     / (period^2);
trend0 = initTrend(data,period);

% level0 = mean(data(1:period));
level0 = mean(data(1:period)) - (period/2)*trend0; % STATA default

season0 = initSeason(data,period);

levelES = zeros(len,1); trendES = zeros(len,1); seasonES = zeros(len,1);
levelES(1) = level0; trendES(1) = trend0; % seasonES(1:freq) = season0;

forecastsIn = zeros(len,1);
forecastsIn(1) = level0 + trend0 + season0(1);
seasonES(1) = gamma * (data(1)-levelES(1)) + (1-gamma) * season0(1);

forecastsOut = zeros(numForecasts,1);
for i=1:(len-1+numForecasts)
    if i >= len
        forecastsOut(i-len+1) = levelES(end) + (i-len+1)*trendES(end) +...
            seasonES(end-period+rem((i+1),period));
    else
        if i < period
            forecastsIn(i+1)...
                = levelES(i) + trendES(i) + season0(i+1);
%             forecastsIn(i+1)...
%                 = (levelES(i) + trendES(i)) * season0(i+1);
            levelES(i+1) = alpha * (data(i+1)-season0(i+1)) +...
                (1-alpha) * (levelES(i)+trendES(i));
%             levelES(i+1) = alpha * (data(i+1)/season0(i+1)) +...
%                 (1-alpha) * (levelES(i)+trendES(i));
            trendES(i+1) = beta * (levelES(i+1)-levelES(i)) +...
                (1-beta) * trendES(i);
            seasonES(i+1) = gamma * (data(i+1)-levelES(i+1)) +...
                (1-gamma) * season0(i+1);
%             seasonES(i+1) = gamma * (data(i+1)/levelES(i+1)) +...
%                 (1-gamma) * season0(i+1);
        else
            forecastsIn(i+1)...
                = levelES(i) + trendES(i) + seasonES(i-period+1);
%             forecastsIn(i+1)...
%                 = (levelES(i) + trendES(i)) * seasonES(i-period+1);
            levelES(i+1) = alpha...
                * (data(i+1)-seasonES(i-period+1))...
                    + (1-alpha) * (levelES(i)+trendES(i));
%             levelES(i+1) = alpha...
%                 * (data(i+1)/seasonES(i-period+1))...
%                     + (1-alpha) * (levelES(i)+trendES(i));
            trendES(i+1) = beta * (levelES(i+1)-levelES(i)) +...
                (1-beta) * trendES(i);
            seasonES(i+1) = gamma * (data(i+1)-levelES(i+1)) +...
                (1-gamma) * seasonES(i-period+1);
%             seasonES(i+1) = gamma * (data(i+1)/levelES(i+1)) +...
%                 (1-gamma) * seasonES(i-period+1);
        end
    end
end


function trend0 = initTrend(data,period) % STATA dafault
len = length(data);

numCycles = floor(len / period);
halfNumCycles = floor(numCycles / 2);

firstAvg = mean(data(1:period));

if numCycles >= 4
    halfAvg =...
        mean(data(((halfNumCycles-1)*period+1):halfNumCycles*period));
    trend0 = (halfAvg-firstAvg) / ((halfNumCycles-1)*period);
elseif and((numCycles >= 2),(numCycles < 4))
    lastAvg =...
        mean(data(((numCycles-1)*period+1):numCycles*period));
    trend0 = (lastAvg-firstAvg) / ((numCycles-1)*period);
elseif and((len > period),(numCycles < 2))
    lastAvg =...
        mean(data((end-period+1):end));
    trend0 = (lastAvg-firstAvg) / (len-period);
else
    trend0 = (data(end)-data(1)) / len;
end

function season0 = initSeason(data,period)
len = length(data);

numCycles = floor(len / period);
cycleAvg = zeros(numCycles,1);
for i=1:numCycles
    cycleAvg(i) = mean(data((period*(i-1)+1):(period*i)));
end

season0 = zeros(period,1);
for i=1:period
    seasonSum = 0;
    for j=1:numCycles
        seasonSum = seasonSum + data(period*(j-1)+i) - cycleAvg(j);
    end
    season0(i) = seasonSum / numCycles;
end


function [forecastsIn,forecastsOut,mse,aic,bic,hqc]...
    = estimateARIMA(data,p,d,q,numForecasts)

if p == 0
    malags = 1:q;
    Mdl = arima('D',d,'MALags',malags);
elseif q == 0
    arlags = 1:p;
    Mdl = arima('D',d,'ARLags',arlags);
else
    arlags = 1:p;
    malags = 1:q;
    Mdl = arima('D',d,'ARLags',arlags,'MALags',malags);
end

% % If heteroscedasticity assumed
% CondVarMdl = garch(1,1);
% Mdl.Variance = CondVarMdl;

[EstMdl,~,logL,~]...
    = estimate(Mdl,data,'Display','off');

resid = infer(EstMdl,data);

forecastsIn = data - resid;
% forecastsIn = max(0, forecastsIn);

[forecastsOut,mse] = forecast(EstMdl,numForecasts,data);

numParam = p + q + 1;
[~,~,ic] = aicbic(logL,numParam,length(data));

aic = ic.aic;
bic = ic.bic;
hqc = ic.hqc;


function [forecastsIn,forecastsOut,mse,aic,bic,hqc]...
    = estimateSARIMA(data,p,d,q,P,D,Q,period,numForecasts)

if (p == 0) && (q == 0)

        if (P == 0)
            smalags = 1:Q;
            % Following MATLAB notation
            Mdl = arima('D',d,'SMALags',period*smalags);

        elseif (Q == 0)
            sarlags = 1:P;
            % Following MATLAB notation
            Mdl = arima('D',d,'SARLags',period*sarlags);
        else
            sarlags = 1:P;
            smalags = 1:Q;
            % Following MATLAB notation
            Mdl = arima('D',d,'SARLags',period*sarlags,'SMALags',period*smalags);
        end

        if D == 1
            Mdl.Seasonality = period;
        end

    % % If heteroscedasticity assumed
    % CondVarMdl = garch(1,1);
    % Mdl.Variance = CondVarMdl;

else
    if (p == 0)
        malags = 1:q;
        
        if (P ~= 0) || (Q ~= 0)

            if (P == 0)
                smalags = 1:Q;
                % Following MATLAB notation
                Mdl = arima('D',d,'MALags',malags,'SMALags',period*smalags);

            elseif (Q == 0)
                sarlags = 1:P;
                % Following MATLAB notation
                Mdl = arima('D',d,'MALags',malags,'SARLags',period*sarlags);
            else
                sarlags = 1:P;
                smalags = 1:Q;
                % Following MATLAB notation
                Mdl = arima('D',d,'MALags',malags,'SARLags',period*sarlags,'SMALags',period*smalags);
            end

            if D == 1
                Mdl.Seasonality = period;
            end

        else
            Mdl = arima('D',d,'MALag',malags);
        end

        % % If heteroscedasticity assumed
        % CondVarMdl = garch(1,1);
        % Mdl.Variance = CondVarMdl;

    elseif (q == 0)
        arlags = 1:p;

        if (P ~= 0) || (Q ~= 0)

            if (P == 0)
                smalags = 1:Q;
                % Following MATLAB notation
                Mdl = arima('ARLags',arlags,'D',d,'SMALags',period*smalags);

            elseif (Q == 0)
                sarlags = 1:P;
                % Following MATLAB notation
                Mdl = arima('ARLags',arlags,'D',d,'SARLags',period*sarlags);
            else
                sarlags = 1:P;
                smalags = 1:Q;
                % Following MATLAB notation
                Mdl = arima('ARLags',arlags,'D',d,'SARLags',period*sarlags,'SMALags',period*smalags);
            end

            if D == 1
                Mdl.Seasonality = period;
            end

        else
            Mdl = arima('ARLags',arlags,'D',d);
        end

        % % If heteroscedasticity assumed
        % CondVarMdl = garch(1,1);
        % Mdl.Variance = CondVarMdl;

    else
        arlags = 1:p;
        malags = 1:q;
        
        if (P ~= 0) || (Q ~= 0)

            if (P == 0)
                smalags = 1:Q;
                % Following MATLAB notation
                Mdl = arima('ARLags',arlags,'D',d,'MALags',malags,'SMALags',period*smalags);

            elseif (Q == 0)
                sarlags = 1:P;
                % Following MATLAB notation
                Mdl = arima('ARLags',arlags,'D',d,'MALags',malags,'SARLags',period*sarlags);
            else
                sarlags = 1:P;
                smalags = 1:Q;
                % Following MATLAB notation
                Mdl = arima('ARLags',arlags,'D',d,'MALags',malags,'SARLags',period*sarlags,'SMALags',period*smalags);
            end

            if D == 1
                Mdl.Seasonality = period;
            end

        else
            Mdl = arima('ARLags',arlags,'D',d,'MALags',malags);
        end

        % % If heteroscedasticity assumed
        CondVarMdl = garch(1,1);
        Mdl.Variance = CondVarMdl;
    end
end

% % If heteroscedasticity assumed
% CondVarMdl = garch(1,1);
% Mdl.Variance = CondVarMdl;

[EstMdl,~,logL,~]...
    = estimate(Mdl,data,'Display','off');

resid = infer(EstMdl,data);

forecastsIn = data - resid;
% forecastsIn = max(0, forecastsIn);

[forecastsOut,mse] = forecast(EstMdl,numForecasts,data);

numParam = p + q + P + Q + 1;
[~,~,ic] = aicbic(logL,numParam,length(data));

aic = ic.aic;
bic = ic.bic;
hqc = ic.hqc;


function [forecastsIn,forecastsOut,level,trend,season,modelName]...
    = forecastFuture(History,modelName,parameters,period,numForecasts)

p = parameters(1);
d = parameters(2);
q = parameters(3);
P = parameters(4);
D = parameters(5);
Q = parameters(6);

[alpha,beta,gamma] = parameterHW(History,period);
[forecastsIn,forecastsOut,level,trend,season]...
    = estimateHW(History,alpha,beta,gamma,period,numForecasts);

if modelName == "HW"

elseif modelName == "ARIMA"
    Err = 1;

    for i = 1:100
        if Err > 0
            try
                [forecastsIn,forecastsOut,~,~,~,~]...
                    = estimateARIMA(History,p,d,q,numForecasts);
                Err = 0;
            catch
                Err = Err + 1;
            end
        end
    end

    if Err > 0
        modelName = "HW";
    end
    
elseif modelName == "SARIMA"
    Err = 1;

    for i = 1:100
        if Err > 0
            try
                [forecastsIn,forecastsOut,~,~,~,~]...
                    = estimateSARIMA(History,p,d,q,P,D,Q,period,numForecasts);
                Err = 0;
            catch
                Err = Err + 1;
            end
        else
            continue
        end
    end

    if Err > 0
        modelName = "HW";
    end

end

forecastsOut = max(0,round(forecastsOut));


%-----------------------------------------------------------------
% Check value of 'dbName' parameter
function [OK,result] = dbNameCheck(dbName)

    if ~or(ischar(dbName),isstring(dbName))
        
        s = struct;
        s.status = "error";
        s.message = "'dbName' must be a char or string type.";

        result = jsonencode(s);

%         fid = fopen('result.json','w');
%         fprintf(fid,'%s',result);
%         fclose(fid);
        
%         error("'dbName' must be a char or string type.");
        
        OK = false;
        
    else
        result = 1;
        
        OK = true;
    end

% Check value of 'dbID' parameter
function [OK,result] = dbIDCheck(dbID)

    if ~or(ischar(dbID),isstring(dbID))
        
        s = struct;
        s.status = "error";
        s.message = "'dbID' must be a char or string type.";

        result = jsonencode(s);

%         fid = fopen('result.json','w');
%         fprintf(fid,'%s',result);
%         fclose(fid);
        
%         error("'dbID' must be a char or string type.");
        
        OK = false;
        
    else
        result = 1;
        
        OK = true;
    end

% Check value of 'dbPW' parameter
function [OK,result] = dbPWCheck(dbPW)

    if ~or(ischar(dbPW),isstring(dbPW))
        
        s = struct;
        s.status = "error";
        s.message = "'dbPW' must be a char or string type.";

        result = jsonencode(s);

%         fid = fopen('result.json','w');
%         fprintf(fid,'%s',result);
%         fclose(fid);
        
%         error("'dbPW' must be a char or string type.");
        
        OK = false;
        
    else
        result = 1;
        
        OK = true;
    end

% Check value of 'VER_ID' parameter
function [OK,result] = VER_IDCheck(VER_ID)

    if ~or(ischar(VER_ID),isstring(VER_ID))
        
        s = struct;
        s.status = "error";
        s.message = "'VER_ID' must be a char or string type.";

        result = jsonencode(s);

%         fid = fopen('result.json','w');
%         fprintf(fid,'%s',result);
%         fclose(fid);
        
%         error("'VER_ID' must be a char or string type.");
        
        OK = false;
        
    else
        result = 1;
        
        OK = true;
    end

% Check value of 'period' and 'fcst' parameters
function [OK,result] = periodCheck(period)

    if ~or(ischar(period),isstring(period))
        
        s = struct;
        s.status = "error";
        s.message = "'periodStart', 'periodEnd', 'fcstStart', and 'fcstEnd' must be char or string type with its length 6 (yyyymm).";

        result = jsonencode(s);

%         fid = fopen('result.json','w');
%         fprintf(fid,'%s',result);
%         fclose(fid);
        
%         error("'DSTRB_CHNL' must be a char or string type.");
        
        OK = false;
        
    elseif strlength(period) ~= 6
        
        s = struct;
        s.status = "error";
        s.message = "'periodStart', 'periodEnd', 'fcstStart', and 'fcstEnd' must be char or string type with its length 6 (yyyymm).";

        result = jsonencode(s);
        
        OK = false;
        
    else
        result = 1;
        
        OK = true;
    end

% Check value of 'DSTRB_CHNL' parameter
function [OK,result] = DSTRB_CHNLCheck(DSTRB_CHNL)

    if ~or(ischar(DSTRB_CHNL),isstring(DSTRB_CHNL))
        
        s = struct;
        s.status = "error";
        s.message = "'DSTRB_CHNL' must be a char or string type.";

        result = jsonencode(s);

%         fid = fopen('result.json','w');
%         fprintf(fid,'%s',result);
%         fclose(fid);
        
%         error("'DSTRB_CHNL' must be a char or string type.");
        
        OK = false;
        
    else
        result = 1;
        
        OK = true;
    end

% Check value of 'SALES_GRP' parameter
function [OK,result] = SALES_GRPCheck(SALES_GRP)

    if ~or(ischar(SALES_GRP),isstring(SALES_GRP))
        
        s = struct;
        s.status = "error";
        s.message = "'SALES_GRP' must be a char or string type.";

        result = jsonencode(s);

%         fid = fopen('result.json','w');
%         fprintf(fid,'%s',result);
%         fclose(fid);
        
%         error("'SALES_GRP' must be a char or string type.");
        
        OK = false;
        
    else
        result = 1;
        
        OK = true;
    end

% Check value of 'parameterInit' parameter
function [OK,result] = parameterInitCheck(parameterInit)

    if ~or(islogical(parameterInit),ismember(parameterInit,[0,1]))
        
        s = struct;
        s.status = "error";
        s.message = "'parameterInit' must be a logical type (true/false or 1/0).";

        result = jsonencode(s);

%         fid = fopen('result.json','w');
%         fprintf(fid,'%s',result);
%         fclose(fid);
        
%         error("'parameterInit' must be a logical type (true/false or 1/0).");
        
        OK = false;
        
    else
        result = 1;
        
        OK = true;
    end

% Check value of 'maxLag' parameter
function [OK,result] = maxLagCheck(maxLag)

    if ~isnumeric(maxLag)
%         error("'maxLag' must be a nonnegative integer.");

        s = struct;
        s.status = "error";
        s.message = "'maxLag' and 'maxLagS' must be nonnegative integers.";

        result = jsonencode(s);
        
        OK = false;
           
    elseif maxLag < 0
%         error("'maxLag' must be a nonnegative integer.");

        s = struct;
        s.status = "error";
        s.message = "'maxLag' and 'maxLagS' must be nonnegative integers.";

        result = jsonencode(s);
        
        OK = false;
        
    else
        result = 1;
        
        OK = true;
    
    end
    
% Check value of 'numTest' parameter
function [OK,result] = numTestCheck(numTest)

    if ~isnumeric(numTest)
%         error("'numTest' must be a nonnegative integer.");

        s = struct;
        s.status = "error";
        s.message = "'numTest' and 'numForecasts' must be nonnegative integers.";

        result = jsonencode(s);
        
        OK = false;
            
    elseif numTest < 0
%         error("'numTest' must be a nonnegative integer.");

        s = struct;
        s.status = "error";
        s.message = "'numTest' and 'numForecasts' must be nonnegative integers.";

        result = jsonencode(s);
        
        OK = false;
        
    else
        result = 1;
        
        OK = true;
    
    end
