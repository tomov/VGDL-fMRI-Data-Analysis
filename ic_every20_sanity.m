% copy of tc_ons_sanity.m

%load mat/ic_ons20.mat
load mat/ic_ons_every20.mat


for r = 1:6
    r

    % interaction changes found when spriteInducting every 20 frames but not normally
    tc_d12 = tc_ons20{r}(find(~ismember(tc_ons20{r}, tc_ons1{r})));
    tc_d12

    ic_d12 = ic_ons20{r}(find(~ismember(ic_ons20{r}, ic_ons1{r})));
    ic_d12

    tec_d12 = tec_ons20{r}(find(~ismember(tec_ons20{r}, tec_ons1{r})));
    tec_d12

    % other way round
    tc_d21_ix = find(~ismember(tc_ons1{r}, tc_ons20{r}));
    tc_d21 = tc_ons1{r}(tc_d21_ix);
    tc_d21
    tc_ons_game_names1{r}(tc_d21_ix)

    ic_d21_ix = find(~ismember(ic_ons1{r}, ic_ons20{r}));
    ic_d21 = ic_ons1{r}(ic_d21_ix);
    ic_d21
    ic_ons_game_names1{r}(ic_d21_ix)

    tec_d21_ix = find(~ismember(tec_ons1{r}, tec_ons20{r}));
    tec_d21 = tec_ons1{r}(tec_d21_ix);
    tec_d21
    tec_ons_game_names1{r}(tec_d21_ix)

end
