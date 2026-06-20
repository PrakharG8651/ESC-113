classdef FluidFlowSimulatorApp < handle

    properties (Access = public)
        UIFigure              matlab.ui.Figure
        LeftPanel             matlab.ui.container.Panel
        ParamPanel            matlab.ui.container.Panel
        NxLabel               matlab.ui.control.Label
        NxEditField           matlab.ui.control.NumericEditField
        NyLabel               matlab.ui.control.Label
        NyEditField           matlab.ui.control.NumericEditField
        DtLabel               matlab.ui.control.Label
        DtEditField           matlab.ui.control.NumericEditField
        MethodPanel           matlab.ui.container.Panel
        MethodGroup           matlab.ui.container.ButtonGroup
        LURadio               matlab.ui.control.RadioButton
        GSRadio               matlab.ui.control.RadioButton
        RunButton             matlab.ui.control.Button
        StatusLabel           matlab.ui.control.Label
        RightPanel            matlab.ui.container.Panel
        PressureAxes          matlab.ui.control.UIAxes
        VelocityAxes          matlab.ui.control.UIAxes
    end

    methods (Access = private)

        function runSimulation(app)
            Nx = app.NxEditField.Value; 
            Ny = app.NyEditField.Value;
            dt = app.DtEditField.Value;
            
            app.StatusLabel.Text = 'Status: Computing...';
            app.StatusLabel.FontColor = [0.9 0.5 0.0];
            drawnow;

            % Grid Setup
            dx = 1 / (Nx - 1);
            dy = 1 / (Ny - 1);
            [X, Y] = meshgrid(linspace(0, 1, Nx), linspace(0, 1, Ny));
            N = Nx * Ny; 

            % Preallocate sparse system
            A = sparse(N, N);
            b = zeros(N, 1);

            % Stencil Assembly
            for j = 1:Ny
                for i = 1:Nx
                    row = i + (j-1)*Nx;
                    if i == 1 || i == Nx || j == 1 || j == Ny
                        A(row, row) = 1;
                        b(row) = 0;
                    else
                        A(row, row) = -2/(dx^2) - 2/(dy^2);
                        A(row, row-1) = 1/(dx^2); 
                        A(row, row+1) = 1/(dx^2); 
                        A(row, row-Nx) = 1/(dy^2); 
                        A(row, row+Nx) = 1/(dy^2); 
                        b(row) = -2 * pi^2 * sin(pi * X(j,i)) * sin(pi * Y(j,i));
                    end
                end
            end

            % Solvers
            if app.LURadio.Value
                [L, U, P_perm] = lu(A); 
                P_vector = U \ (L \ (P_perm * b));
            else
                P_vector = zeros(N, 1); 
                maxIter = 1000;
                tol = 1e-5;
                D = diag(diag(A));
                L_mat = tril(A, -1);
                U_mat = triu(A, 1);
                
                for iter = 1:maxIter
                    P_old = P_vector;
                    P_vector = (D + L_mat) \ (b - U_mat * P_old);
                    if norm(P_vector - P_old, inf) < tol
                        break;
                    end
                end
            end

            P = reshape(P_vector, [Ny, Nx]);
            [dPdx, dPdy] = gradient(P, dx, dy);
            U_vel = -dt * dPdx;
            V_vel = -dt * dPdy;

            % Render Pressure
            cla(app.PressureAxes);
            contourf(app.PressureAxes, X, Y, P, 20, 'LineColor', 'none');
            colorbar(app.PressureAxes);
            colormap(app.PressureAxes, 'turbo');
            title(app.PressureAxes, 'Pressure Field (P)');
            xlabel(app.PressureAxes, 'X'); ylabel(app.PressureAxes, 'Y');

            % Render Velocity
            cla(app.VelocityAxes);
            quiver(app.VelocityAxes, X, Y, U_vel, V_vel, 1.5, 'r', 'LineWidth', 1);
            title(app.VelocityAxes, 'Velocity Field Vectors (U, V)');
            xlabel(app.VelocityAxes, 'X'); ylabel(app.VelocityAxes, 'Y');
            app.VelocityAxes.XLim = [0 1]; app.VelocityAxes.YLim = [0 1];

            app.StatusLabel.Text = 'Status: Ready';
            app.StatusLabel.FontColor = [0.1 0.6 0.1];
        end

        function setupUI(app)
            % Clean configuration window frame
            app.UIFigure = uifigure('Name', 'Fluid Flow Simulator Pro', 'Position', [100, 100, 950, 520]);
            app.UIFigure.Color = [0.94 0.95 0.97];

            % Control Column
            app.LeftPanel = uipanel(app.UIFigure, 'Title', 'Control Panel', 'Position', [15, 15, 260, 490], 'FontWeight', 'bold', 'BackgroundColor', 'white');

            % Parameters Group
            app.ParamPanel = uipanel(app.LeftPanel, 'Title', '1. Discretization Settings', 'Position', [10, 270, 238, 190], 'BackgroundColor', 'white');
            app.NxLabel = uilabel(app.ParamPanel, 'Position', [15, 130, 100, 22], 'Text', 'Nodes X (Nx):');
            app.NxEditField = uieditfield(app.ParamPanel, 'numeric', 'Position', [130, 130, 80, 22], 'Value', 20, 'LowerLimit', 5, 'UpperLimit', 100, 'RoundFractionalValues', 'on');
            app.NyLabel = uilabel(app.ParamPanel, 'Position', [15, 85, 100, 22], 'Text', 'Nodes Y (Ny):');
            app.NyEditField = uieditfield(app.ParamPanel, 'numeric', 'Position', [130, 85, 80, 22], 'Value', 20, 'LowerLimit', 5, 'UpperLimit', 100, 'RoundFractionalValues', 'on');
            app.DtLabel = uilabel(app.ParamPanel, 'Position', [15, 40, 100, 22], 'Text', 'Time Step (dt):');
            app.DtEditField = uieditfield(app.ParamPanel, 'numeric', 'Position', [130, 40, 80, 22], 'Value', 0.1, 'LowerLimit', 0.001, 'UpperLimit', 10);

            % Solver Group
            app.MethodPanel = uipanel(app.LeftPanel, 'Title', '2. Numerical Engine', 'Position', [10, 130, 238, 120], 'BackgroundColor', 'white');
            app.MethodGroup = uibuttongroup(app.MethodPanel, 'BorderType', 'none', 'Position', [10, 10, 210, 80], 'BackgroundColor', 'white');
            app.LURadio = uiradiobutton(app.MethodGroup, 'Text', 'LU Decomposition', 'Position', [10, 45, 180, 22], 'Value', true);
            app.GSRadio = uiradiobutton(app.MethodGroup, 'Text', 'Gauss-Seidel Iterative', 'Position', [10, 15, 180, 22]);

            % Run Button
            app.RunButton = uibutton(app.LeftPanel, 'push', 'Position', [10, 60, 238, 45], 'Text', 'RUN SIMULATION', ...
                'FontWeight', 'bold', 'FontSize', 13, 'BackgroundColor', [0.0 0.45 0.74], 'FontColor', 'white', ...
                'ButtonPushedFcn', @(btn, event) app.runSimulation());

            % Status Output
            app.StatusLabel = uilabel(app.LeftPanel, 'Position', [15, 20, 230, 22], 'Text', 'Status: Ready', 'FontWeight', 'bold', 'FontColor', [0.1 0.6 0.1]);

            % Viewport Display Grid Panel
            app.RightPanel = uipanel(app.UIFigure, 'Title', 'Simulation Viewports', 'Position', [290, 15, 645, 490], 'FontWeight', 'bold', 'BackgroundColor', 'white');
            app.PressureAxes = uiaxes(app.RightPanel, 'Position', [15, 45, 295, 390]);
            app.VelocityAxes = uiaxes(app.RightPanel, 'Position', [330, 45, 295, 390]);
        end
    end

    methods (Access = public)
        function app = FluidFlowSimulatorApp()
            setupUI(app);
            runSimulation(app); 
        end
    end
end