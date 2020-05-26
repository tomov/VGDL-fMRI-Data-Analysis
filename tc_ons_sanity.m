load mat/tc_ons.mat;

% ver 1: before refactor, full finalTimeStepList history
% ver 2: before refactor, reset finalTimeStepList for each play
% ver 3: refactor, reset finalTimeStepList for each play

% almost all the differences are explained by sprites, except for one in run 4

% =>

for r = 5:6
    r

    tc_d12 = tc_ons2{r}(find(~ismember(tc_ons2{r}, tc_ons3{r})));
    tc_d21 = tc_ons3{r}(find(~ismember(tc_ons3{r}, tc_ons2{r})));
    sc_d12 = sc_ons2{r}(find(~ismember(sc_ons2{r}, sc_ons3{r})));
    sc_d21 = sc_ons3{r}(find(~ismember(sc_ons3{r}, sc_ons2{r})));
    ic_d12 = ic_ons2{r}(find(~ismember(ic_ons2{r}, ic_ons3{r})));
    ic_d21 = ic_ons3{r}(find(~ismember(ic_ons3{r}, ic_ons2{r})));
    tec_d12 = tec_ons2{r}(find(~ismember(tec_ons2{r}, tec_ons3{r})));
    tec_d21 = tec_ons3{r}(find(~ismember(tec_ons3{r}, tec_ons2{r})));

    if length(tc_d12) > 0
        immse(tc_d12, sc_d12)
    end
    if length(tc_d21) > 0
        immse(tc_d21, sc_d21)
    end
end

%{
for r = 5:6
    r

    tc_d23 = tc_ons2{r}(find(~ismember(tc_ons2{r}, tc_ons3{r})));
    tc_d32 = tc_ons3{r}(find(~ismember(tc_ons3{r}, tc_ons2{r})));
    sc_d23 = sc_ons2{r}(find(~ismember(sc_ons2{r}, sc_ons3{r})));
    sc_d32 = sc_ons3{r}(find(~ismember(sc_ons3{r}, sc_ons2{r})));
    ic_d23 = ic_ons2{r}(find(~ismember(ic_ons2{r}, ic_ons3{r})));
    ic_d32 = ic_ons3{r}(find(~ismember(ic_ons3{r}, ic_ons2{r})));
    tec_d23 = tec_ons2{r}(find(~ismember(tec_ons2{r}, tec_ons3{r})));
    tec_d32 = tec_ons3{r}(find(~ismember(tec_ons3{r}, tec_ons2{r})));

    if length(tc_d23) > 0
        immse(tc_d23, sc_d23)
    end
    if length(tc_d32) > 0
        immse(tc_d32, sc_d32)
    end
end
%}
