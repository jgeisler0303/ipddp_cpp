function ipddp_problem_cpp(problem_name, problem_func)
% problem_name= 'Car';
[funcs, fp, ~, sym_vars] = problem_func();

nx= length(sym_vars.x);
nu= length(sym_vars.u);
nc= length(sym_vars.c);

[const, nonconst]= detect_const(funcs, nx+nu);

ll= strings(0);

ll(1, 1)= "#include <Eigen/Dense>";
ll(2, 1)= "";
ll(end+1)=  string(sprintf("constexpr int nx= %d;", nx));
ll(end+1)=  string(sprintf("constexpr int nu= %d;", nu));
ll(end+1)=  string(sprintf("constexpr int nc= %d;", nc));
ll(end+1)=  string(sprintf("constexpr int has_fxx= %d;", isfield(sym_vars, 'fxx')));
ll(end+1)=  string(sprintf("constexpr int has_cxx= %d;", isfield(sym_vars, 'cxx')));
ll(end+1)=  string(sprintf("constexpr int n_hor= %d;", fp.horizon+1));
ll(end+1)= "";
ll(end+1)=  "typedef Eigen::Matrix<double, nx, 1> VecX;";
ll(end+1)=  "typedef Eigen::Matrix<double, nu, 1> VecU;";
ll(end+1)=  "typedef Eigen::Matrix<double, nc, 1> VecC;";
ll(end+1)=  "typedef Eigen::Matrix<double, nx, nx> MatXX;";
ll(end+1)=  "typedef Eigen::Matrix<double, nu, nu> MatUU;";
ll(end+1)=  "typedef Eigen::Matrix<double, nx, nu> MatXU;";
ll(end+1)=  "typedef Eigen::Matrix<double, nu, nx> MatUX;";
ll(end+1)=  "typedef Eigen::Matrix<double, nc, nx> MatCX;";
ll(end+1)=  "typedef Eigen::Matrix<double, nc, nu> MatCU;";
ll(end+1)=  "typedef Eigen::Matrix<double, nc, nc> MatCC;";
ll(end+1)=  "typedef std::array<MatXX, nx> TensXXX;";
ll(end+1)=  "typedef std::array<MatXU, nx> TensXXU;";
ll(end+1)=  "typedef std::array<MatUU, nx> TensXUU;";
ll(end+1)=  "typedef std::array<MatXX, nc> TensCXX;";
ll(end+1)=  "typedef std::array<MatXU, nc> TensCXU;";
ll(end+1)=  "typedef std::array<MatUU, nc> TensCUU;";
ll(end+1)= "";
ll(end+1)= "class Problem"+problem_name+" {";
ll(end+1)= "public:";
ll(end+1)= "    Problem"+problem_name+"() {";
ll(end+1)= "        x.setZero();";
ll(end+1)= "        u.setZero();";
ll(end+1)= "        c.setZero();";
ll(end+1)= "        y.setZero();";
ll(end+1)= "        s.setZero();";
% ll(end+1)= "        mu(4).setZero();";
ll(end+1)= "        px.setZero();"; 
ll(end+1)= "        pxx.setZero();";
ll(end+1)= "        fx.setZero();";
ll(end+1)= "        fu.setZero();";
if isfield(sym_vars, 'fxx')
    ll(end+1)= "        for(int i= 0; i<nc; ++i) {";
    ll(end+1)= "            fxx[i].setZero();";
    ll(end+1)= "            fuu[i].setZero();";
    ll(end+1)= "            fxu[i].setZero();";
    ll(end+1)= "        }";
end
ll(end+1)= "        cx.setZero();";
ll(end+1)= "        cu.setZero();";
if isfield(sym_vars, 'cxx')
    ll(end+1)= "        for(int i= 0; i<nc; ++i) {";
    ll(end+1)= "            cxx[i].setZero();";
    ll(end+1)= "            cuu[i].setZero();";
    ll(end+1)= "            cxu[i].setZero();";
    ll(end+1)= "        }";
end
ll(end+1)= "        qx.setZero();";
ll(end+1)= "        qu.setZero();";
ll(end+1)= "        qxx.setZero();";
ll(end+1)= "        qxu.setZero();";
ll(end+1)= "        quu.setZero();";
ll(end+1)= "        ku.setZero();";
ll(end+1)= "        ky.setZero();";
ll(end+1)= "        ks.setZero();";
ll(end+1)= "        Ku.setZero();";
ll(end+1)= "        Ky.setZero();";
ll(end+1)= "        Ks.setZero();";
ll(end+1)= "";
ll= [ll; set_const(const)];
ll(end+1)= "    }";
ll(end+1)= "";

ll(end+1)= "    VecX step() {";
ll(end+1)= "        VecX f;";
tmp_file= [tempname '.m'];
matlabFunction(sym_vars.f, 'Vars', {sym_vars.x, sym_vars.u}, 'Outputs', {'f'}, 'File', tmp_file);
ll_next= parse_calculations(tmp_file, sym_vars, nonconst);
delete(tmp_file)
ll= [ll; ll_next];
ll(end+1)= "        return f;";
ll(end+1)= "    }";
ll(end+1)= "";

ll(end+1)= "    void calc_c() {";
tmp_file= [tempname '.m'];
matlabFunction(sym_vars.c, 'Vars', {sym_vars.x, sym_vars.u}, 'Outputs', {'c'}, 'File', tmp_file);
ll_next= parse_calculations(tmp_file, sym_vars, nonconst);
delete(tmp_file)
ll= [ll; ll_next];
ll(end+1)= "    }";
ll(end+1)= "";

ll(end+1)= "    double calc_q() {";
tmp_file= [tempname '.m'];
matlabFunction(sym_vars.q, 'Vars', {sym_vars.x, sym_vars.u}, 'Outputs', {'q'}, 'File', tmp_file);
ll_next= parse_calculations(tmp_file, sym_vars, nonconst);
delete(tmp_file)
ll= [ll; ll_next];
ll(end+1)= "        return q;";
ll(end+1)= "    }";
ll(end+1)= "";

ll(end+1)= "    double calc_p() {";
tmp_file= [tempname '.m'];
matlabFunction(sym_vars.p, 'Vars', {sym_vars.x, sym_vars.u}, 'Outputs', {'p'}, 'File', tmp_file);
ll_next= parse_calculations(tmp_file, sym_vars, nonconst);
delete(tmp_file)
ll= [ll; ll_next];
ll(end+1)= "        return p;";
ll(end+1)= "    }";
ll(end+1)= "";

ll(end+1)= "    void calc_grads() {";
tmp_file= [tempname '.m'];
matlabFunction(sym_vars.fx, sym_vars.fxx, sym_vars.fu, sym_vars.fuu, sym_vars.fxu, sym_vars.cx, sym_vars.cu, sym_vars.qx, sym_vars.qu, sym_vars.qxx, sym_vars.quu, sym_vars.qxu, 'Vars', {sym_vars.x, sym_vars.u}, 'Outputs', {'fx', 'fxx', 'fu', 'fuu', 'fxu', 'cx', 'cu', 'qx', 'qu', 'qxx', 'quu', 'qxu'}, 'File', tmp_file);
ll_grads= parse_calculations(tmp_file, sym_vars, nonconst);
delete(tmp_file)
ll= [ll; ll_grads];
ll(end+1)= "    }";
ll(end+1)= "";

ll(end+1)= "    void calc_final_grads() {";
tmp_file= [tempname '.m'];
matlabFunction(sym_vars.px, sym_vars.pxx, 'Vars', {sym_vars.x, sym_vars.u}, 'Outputs', {'px', 'pxx'}, 'File', tmp_file);
ll_grads= parse_calculations(tmp_file, sym_vars, nonconst);
delete(tmp_file)
ll= [ll; ll_grads];
ll(end+1)= "    }";
ll(end+1)= "";
ll(end+1)= "    double q, p;";
ll(end+1)= "    VecX x;";
ll(end+1)= "    VecU u;";
ll(end+1)= "    VecC c;";
ll(end+1)= "    VecC y;";
ll(end+1)= "    VecC s;";
% ll(end+1)= "    VecC mu;";
ll(end+1)= "    ";
ll(end+1)= "    VecX px;  // Assuming px is the gradient of p with respect to x";
ll(end+1)= "    MatXX pxx; // Assuming pxx is the Hessian of p with respect to x";
ll(end+1)= "    ";
ll(end+1)= "    MatXX fx;";
ll(end+1)= "    MatXU fu;";
ll(end+1)= "";
if isfield(sym_vars, 'fxx')
    ll(end+1)= "    TensXXX fxx;";
    ll(end+1)= "    TensXUU fuu;";
    ll(end+1)= "    TensXXU fxu;";
    ll(end+1)= "";
end
ll(end+1)= "    MatCX cx;";
ll(end+1)= "    MatCU cu;";
if isfield(sym_vars, 'cxx')
    ll(end+1)= "    TensCXX cxx;";
    ll(end+1)= "    TensCUU cuu;";
    ll(end+1)= "    TensCXU cxu;";
    ll(end+1)= "";
end
ll(end+1)= "    ";
ll(end+1)= "    VecX qx; ";
ll(end+1)= "    VecU qu;";
ll(end+1)= "    MatXX qxx;";
ll(end+1)= "    MatXU qxu;";
ll(end+1)= "    MatUU quu;";
ll(end+1)= "";
ll(end+1)= "    VecU ku;";
ll(end+1)= "    VecC ky;";
ll(end+1)= "    VecC ks;";
ll(end+1)= "    MatUX Ku;";
ll(end+1)= "    MatCX Ky;";
ll(end+1)= "    MatCX Ks;";
ll(end+1)= "};";

writelines(ll, ['ipddp_problem_' problem_name '.hpp'])
end

function ll_out= parse_calculations(file_name, sym_vars, nonconst)
indent= "        ";
ll= readlines(file_name, 'WhitespaceRule', 'trim', 'EmptyLineRule', 'skip');

for i= length(ll):-1:1
    if startsWith(ll(i), 'if') || ...
        startsWith(ll(i), 'end') || ...
        startsWith(ll(i), '%') || ...
        startsWith(ll(i), 'function') || ...
        startsWith(ll(i), 'u') || ...
        startsWith(ll(i), 'x')
        
        ll(i)= [];
    end
end

ll_out= strings(0);
for i= 1:length(ll)
    l= ll(i);
    
    if(contains(l, "reshape"))
        lla= assign_matrix(l, nonconst);
        ll_out= [ll_out; lla'];
    elseif(contains(l, "["))
        lla= assign_vector(l, nonconst);
        ll_out= [ll_out; lla'];
    else
        ll_out= [ll_out; l];
    end
end
for i= 1:length(ll_out)
    l= ll_out(i);
    l= strrep(l, '.*', '*');
    l= strrep(l, './', '/');
    l= regexprep(l, '\<pi\>', 'M_PI');

    ind_pows= strfind(l, '.^');
    l_pow= l;
    for ip= length(ind_pows):-1:1
        open_paren= 0;
        for i_base= (ind_pows(ip)-1):-1:1
            if l{1}(i_base)==')'
                open_paren= open_paren + 1;
            elseif l{1}(i_base)=='('
                open_paren= open_paren - 1;
            end

            if open_paren<=0 && any(l{1}(i_base)=='=+-*/( ')
                if l{1}(i_base)=='+' && i_base>1 && lower(l{1}(i_base-1))=='e'
                else
                    break;
                end
            end
        end
        base= l{1}((i_base+1):(ind_pows(ip)-1));

        open_paren= 0;
        for i_expo= (ind_pows(ip)+2):strlength(l)
            if l{1}(i_expo)==')'
                open_paren= open_paren - 1;
            elseif l{1}(i_expo)=='('
                open_paren= open_paren + 1;
            end

            if open_paren<=0 && any(l{1}(i_expo)=='+-*/;) ')
                if l{1}(i_expo)=='+' && lower(l{1}(i_expo-1))=='e'
                else
                    break;
                end
            end
        end
        expo= l{1}((ind_pows(ip)+2):(i_expo-1));
        l_pow= strrep(l_pow, l{1}((i_base+1):(i_expo-1)), sprintf('pow(%s, %s)', base, expo));
    end
    l= l_pow;

    if(startsWith(l, 't'))
        l= "double " + l;
    else
        var_expr= strsplit(l, '=');
        var_expr(1)= strip(var_expr(1));
        if strlength(var_expr(1))>1 && (var_expr{1}(end)=='x' || var_expr{1}(end)=='u')
            l= var_expr(1)+"(0)= "+var_expr(2);
        end
    end
    for ix= 1:length(sym_vars.x)
        l= regexprep(l, ['\<' char(sym_vars.x(ix)) '\>'], num2str(ix-1, 'x(%d)'));
    end
    for iu= 1:length(sym_vars.u)
        l= regexprep(l, ['\<' char(sym_vars.u(iu)) '\>'], num2str(iu-1, 'u(%d)'));
    end

    ll_out(i)= l;
end
ll_out= append(indent, ll_out);
end

function lla= assign_matrix(l, nonconst)
t= regexp(l, '(.*) = reshape\(\[(.*?)\],\[(.*?)\]\);', 'tokens');
elements= strsplit(t{1}{2}, ',');
dims= str2double(strsplit(t{1}{3}, ','));
elements= reshape(elements, dims);
var= t{1}{1};

lla= strings(0);
if length(dims)==2
    for row= 1:dims(1)
        for col= 1:dims(2)
            if nonconst.(var)(row, col)
                lla(end+1)= string(sprintf("%s(%d, %d)= %s;", var, row-1, col-1, elements{row, col}));
            end
        end
    end
else
    for elem= 1:dims(1)
        for row= 1:dims(2)
            for col= 1:dims(3)
                if nonconst.(var)(elem, row, col)
                    lla(end+1)= string(sprintf("%s[%d](%d, %d)= %s;", var, elem-1, row-1, col-1, elements{elem, row, col}));
                end
            end
        end
    end
end
end

function lla= assign_vector(l, nonconst)
t= regexp(l, '(.*) = \[(.*?)\];', 'tokens');
elements= strsplit(t{1}{2}, {';' ','});
var= t{1}{1};

lla= strings(0);
for row= 1:length(elements)
    if nonconst.(var)(row)
        lla(end+1)= string(sprintf("%s(%d)= %s;", var, row-1, elements{row}));
    end
end
end

function [const, nonconst]= detect_const(funcs, nxu)
nans= nan(nxu, 1);
fields= fieldnames(funcs);
for i= 1:length(fields)
    f= funcs.(fields{i});
    if ~isa(f, 'function_handle'), continue, end
    func= f(nans);
    nan_func= isnan(func);
    nonconst.(fields{i})= nan_func;
    func(nan_func)= 0;
    const.(fields{i})= func;
end
end

function ll= set_const(consts)
ll= strings(0);
fields= fieldnames(consts);
for i= 1:length(fields)
    dims= size(consts.(fields{i}));
    if length(dims)==2
        dims_= [1 dims];
    else
        dims_= dims;
    end
    for elem= 1:dims_(1)
        for row= 1:dims_(2)
            for col= 1:dims_(3)
                if length(dims)==2
                    c= consts.(fields{i})(row, col);
                    if c==0, continue, end
                    if dims(2)==1
                        ll(end+1)= string(sprintf('        %s(%d)= %g;', fields{i}, row-1, c));
                    else
                        ll(end+1)= string(sprintf('        %s(%d, %d)= %g;', fields{i}, row-1, col-1, consts.(fields{i})(row, col)));
                    end
                else
                    c= consts.(fields{i})(row, col);
                    if c==0, continue, end
                    ll(end+1)= string(sprintf('        %s[%d](%d)= %g;', fields{i}, elem-1, row-1, col-1, c));
                end
            end
        end
    end
    
end
ll= ll';
end