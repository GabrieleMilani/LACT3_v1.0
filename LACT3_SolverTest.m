%-------------------------------------------------------------------------%
%           LIMIT ANALYSIS CODE FOR TILTING TABLE TESTS (LACT3)           %
%                        Automatic solver test                            %
%-------------------------------------------------------------------------%

function tests = LACT3_SolverTest
    clear all;
    clc;
    warning off;
    tests = functiontests(localfunctions);
end

function testPortal(testCase)
    pathname="sample\";
    filename="Portal.dxf";
    teta=27.30*pi/180;
    c=0;
    phi=30;
    [P, polyL, int_i] = Load_wall_data_2(filename, pathname,[]);
    [polyL, int_i] = Material(polyL, int_i, c, phi, 1, 25);
    [~,ub] = Solver_in_UB(P, polyL, int_i, teta);
    [~,lb] = Solver_in_LB(P, polyL, int_i, teta);
    actSolution = [abs(ub)<1e-3 abs(lb)<1e-3];
    expSolution = [true true];
    verifyEqual(testCase,actSolution,expSolution)
end
% 
function testwall(testCase)
    pathname="sample\";
    filename="wall.dxf";
    teta=16.73*pi/180;
    c=0;
    phi=26;
    [P, polyL, int_i] = Load_wall_data_2(filename, pathname,[]);
    [polyL, int_i] = Material(polyL, int_i, c, phi, 1, 25);
    [~,ub] = Solver_in_UB(P, polyL, int_i, teta);
    [~,lb] = Solver_in_LB(P, polyL, int_i, teta);
    actSolution = [abs(ub)<1e-3 abs(lb)<1e-3];
    expSolution = [true true];
    verifyEqual(testCase,actSolution,expSolution)
end
% 
function testarch_1(testCase)
    pathname="sample\";
    filename="arch_1.dxf";
    teta=17.10*pi/180;
    c=0;
    phi=30;
    [P, polyL, int_i] = Load_wall_data_2(filename, pathname,[]);
    [polyL, int_i] = Material(polyL, int_i, c, phi, 1, 25);
    [~,ub] = Solver_in_UB(P, polyL, int_i, teta);
    [~,lb] = Solver_in_LB(P, polyL, int_i, teta);
    actSolution = [abs(ub)<1e-3 abs(lb)<1e-3];
    expSolution = [true true];
    verifyEqual(testCase,actSolution,expSolution)
end