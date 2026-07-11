pub trait Shape1D {
    // quadrature points
    fn num_quad(&self) -> usize;
    fn w(&self) -> Vec<f64>;
    fn a(&self) -> Vec<f64>;

    // shape functions
    fn num_node(&self) -> usize;
    fn eval(&self, a: f64) -> Vec<f64>;
    fn grad(&self, a: f64) -> Vec<f64>;

    // shape function evaluation
    // access using [q][v]; value for shape function v at quadrature point q
    fn n(&self) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
        // quadrature points
        let num_quad = self.num_quad();
        let quad_a = self.a();

        // shape functions
        let mut quad_n = Vec::with_capacity(num_quad);
        let mut quad_gna = Vec::with_capacity(num_quad);
        for qid in 0..num_quad {
            // values at integration point
            let n = self.eval(quad_a[qid]);
            let dnda = self.grad(quad_a[qid]);

            // store values
            quad_n.push(n);
            quad_gna.push(dnda);
        }

        // return
        (quad_n, quad_gna)
    }

}

pub trait Shape2D {
    // quadrature points
    fn num_quad(&self) -> usize;
    fn w(&self) -> Vec<f64>;
    fn a(&self) -> Vec<f64>;
    fn b(&self) -> Vec<f64>;

    // shape functions
    fn num_node(&self) -> usize;
    fn eval(&self, a: f64, b: f64) -> Vec<f64>;
    fn grad(&self, a: f64, b: f64) -> (Vec<f64>, Vec<f64>);

    // shape function evaluation
    // access using [q][v]; value for shape function v at quadrature point q
    fn n(&self) -> (Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<Vec<f64>>) {
        // quadrature points
        let num_quad = self.num_quad();
        let quad_a = self.a();
        let quad_b = self.b();

        // shape functions
        let mut quad_n = Vec::with_capacity(num_quad);
        let mut quad_gna = Vec::with_capacity(num_quad);
        let mut quad_gnb = Vec::with_capacity(num_quad);
        for qid in 0..num_quad {
            // values at integration point
            let n = self.eval(quad_a[qid], quad_b[qid]);
            let (dnda, dndb) = self.grad(quad_a[qid], quad_b[qid]);

            // store values
            quad_n.push(n);
            quad_gna.push(dnda);
            quad_gnb.push(dndb);
        }

        // return
        (quad_n, quad_gna, quad_gnb)
    }
}
