function projOp(ϵ, f, ∂f)
    eigenValues, eigenVectors = eigen(ϵ)
    delta_lambda = eigenValues[1] - eigenValues[2]

    eigenVectors_1o1 = eigenVectors[:, 1] * eigenVectors[:, 1]'
    eigenVectors_2o2 = eigenVectors[:, 2] * eigenVectors[:, 2]'
    P1 = ∂f_function(eigenValues[1]) * (eigenVectors_1o1 ⊗ eigenVectors_1o1) +
         ∂f_function(eigenValues[2]) * (eigenVectors_2o2 ⊗ eigenVectors_2o2)

    if abs(delta_lambda) < eps()
        return P1
    else
        eigenVectors_1o2 = eigenVectors[:, 1] * eigenVectors[:, 2]'
        eigenVectors_1o2p2o1 = eigenVectors_1o2 + eigenVectors_1o2'

        factor = 0.5 * ((f_function(eigenValues[1]) - f_function(eigenValues[2])) / delta_lambda)
        term = eigenVectors_1o2p2o1 ⊗ eigenVectors_1o2p2o1
        P2 = factor * term

        return P1 + P2
    end
end