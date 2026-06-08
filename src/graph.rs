use crate::error::{WaveError, WaveResult};
use serde::{Deserialize, Serialize};

/// An undirected graph stored as an adjacency list.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Graph {
    adj: Vec<Vec<usize>>,
    n: usize,
}

impl Graph {
    /// Create an empty graph with `n` nodes and no edges.
    pub fn new(n: usize) -> Self {
        Self {
            adj: vec![vec![]; n],
            n,
        }
    }

    /// Number of nodes.
    pub fn n(&self) -> usize {
        self.n
    }

    /// Number of (undirected) edges.
    pub fn m(&self) -> usize {
        self.adj.iter().map(|v| v.len()).sum::<usize>() / 2
    }

    /// Add an undirected edge between `u` and `v`. No-op if already present.
    pub fn add_edge(&mut self, u: usize, v: usize) -> WaveResult<()> {
        if u >= self.n || v >= self.n {
            return Err(WaveError::IndexOutOfBounds {
                index: u.max(v),
                len: self.n,
            });
        }
        if u == v {
            return Ok(());
        }
        if !self.adj[u].contains(&v) {
            self.adj[u].push(v);
            self.adj[v].push(u);
        }
        Ok(())
    }

    /// Neighbours of node `i`.
    pub fn neighbors(&self, i: usize) -> &[usize] {
        &self.adj[i]
    }

    /// Degree of node `i`.
    pub fn degree(&self, i: usize) -> usize {
        self.adj[i].len()
    }

    /// Build the combinatorial Laplacian `L = D - A` as a dense matrix.
    pub fn laplacian(&self) -> Vec<Vec<f64>> {
        let n = self.n;
        let mut l = vec![vec![0.0; n]; n];
        for (i, row) in l.iter_mut().enumerate() {
            row[i] = self.adj[i].len() as f64;
            for &j in &self.adj[i] {
                row[j] -= 1.0;
            }
        }
        l
    }

    /// Compute all eigenvalues and eigenvectors of the Laplacian using the
    /// Jacobi eigenvalue algorithm for symmetric matrices.
    fn full_eigen(&self) -> (Vec<f64>, Vec<Vec<f64>>) {
        let n = self.n;
        if n == 0 {
            return (vec![], vec![]);
        }
        let mut a = self.laplacian();
        let mut v = identity(n);

        for _ in 0..100 * n * n {
            // Find the largest off-diagonal element
            let mut max_val = 0.0_f64;
            let mut p = 0;
            let mut q = 1;
            for (i, row) in a.iter().enumerate() {
                for (j, val) in row.iter().enumerate().skip(i + 1) {
                    if val.abs() > max_val {
                        max_val = val.abs();
                        p = i;
                        q = j;
                    }
                }
            }
            if max_val < 1e-14 {
                break;
            }

            // Compute rotation angle
            let app = a[p][p];
            let aqq = a[q][q];
            let apq = a[p][q];
            let theta = if (app - aqq).abs() < 1e-30 {
                std::f64::consts::FRAC_PI_4
            } else {
                0.5 * (2.0 * apq / (app - aqq)).atan()
            };
            let c = theta.cos();
            let s = theta.sin();

            // Apply rotation
            for i in 0..n {
                if i != p && i != q {
                    let aip = a[i][p];
                    let aiq = a[i][q];
                    a[i][p] = c * aip + s * aiq;
                    a[p][i] = a[i][p];
                    a[i][q] = -s * aip + c * aiq;
                    a[q][i] = a[i][q];
                }
                let vip = v[i][p];
                let viq = v[i][q];
                v[i][p] = c * vip + s * viq;
                v[i][q] = -s * vip + c * viq;
            }
            let new_pp = c * c * app + 2.0 * s * c * apq + s * s * aqq;
            let new_qq = s * s * app - 2.0 * s * c * apq + c * c * aqq;
            a[p][p] = new_pp;
            a[q][q] = new_qq;
            a[p][q] = 0.0;
            a[q][p] = 0.0;
        }

        // Extract eigenvalues and eigenvectors
        let mut eigenvalues: Vec<(f64, usize)> = (0..n).map(|i| (a[i][i], i)).collect();
        eigenvalues.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        let vals: Vec<f64> = eigenvalues.iter().map(|(v, _)| *v).collect();
        // Eigenvectors are columns of v, reorder
        let perm: Vec<usize> = eigenvalues.iter().map(|(_, i)| *i).collect();
        let vecs: Vec<Vec<f64>> = (0..n)
            .map(|i| (0..n).map(|j| v[i][perm[j]]).collect())
            .collect();

        // vecs[j] is the j-th eigenvector (as column)
        let eigvecs: Vec<Vec<f64>> = (0..n)
            .map(|k| (0..n).map(|i| vecs[i][k]).collect())
            .collect();

        (vals, eigvecs)
    }

    /// Compute the second-smallest eigenvalue `λ₂` (algebraic connectivity / Fiedler value).
    pub fn fiedler(&self) -> WaveResult<(f64, Vec<f64>)> {
        if self.n == 0 {
            return Err(WaveError::EmptyGraph);
        }
        if self.n == 1 {
            return Ok((0.0, vec![1.0]));
        }

        let (vals, vecs) = self.full_eigen();

        // λ₂ is the second eigenvalue (index 1)
        let lambda2 = vals[1];
        if lambda2 < 1e-10 {
            return Err(WaveError::Disconnected);
        }

        let mut fiedler_vec = vecs[1].clone();
        // Normalize
        let norm = fiedler_vec.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm > 1e-15 {
            for x in fiedler_vec.iter_mut() {
                *x /= norm;
            }
        }

        Ok((lambda2, fiedler_vec))
    }

    /// Compute the `k` smallest non-zero eigenvalues and eigenvectors.
    pub fn eigenmodes(&self, k: usize) -> WaveResult<Vec<(f64, Vec<f64>)>> {
        if self.n == 0 {
            return Err(WaveError::EmptyGraph);
        }

        let (vals, vecs) = self.full_eigen();

        // Skip eigenvalue 0 (index 0), take next k
        let modes: Vec<(f64, Vec<f64>)> = (1..=(k.min(self.n - 1)))
            .map(|i| {
                let mut v = vecs[i].clone();
                let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
                if norm > 1e-15 {
                    for x in v.iter_mut() {
                        *x /= norm;
                    }
                }
                (vals[i], v)
            })
            .collect();

        Ok(modes)
    }

    /// Check if the graph is connected using BFS.
    pub fn is_connected(&self) -> bool {
        if self.n == 0 {
            return true;
        }
        let mut visited = vec![false; self.n];
        let mut stack = vec![0];
        visited[0] = true;
        let mut count = 1;
        while let Some(u) = stack.pop() {
            for &v in &self.adj[u] {
                if !visited[v] {
                    visited[v] = true;
                    stack.push(v);
                    count += 1;
                }
            }
        }
        count == self.n
    }
}

fn identity(n: usize) -> Vec<Vec<f64>> {
    let mut m = vec![vec![0.0; n]; n];
    for (i, row) in m.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    m
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_graph() {
        let g = Graph::new(0);
        assert_eq!(g.n(), 0);
        assert_eq!(g.m(), 0);
    }

    #[test]
    fn test_single_node() {
        let g = Graph::new(1);
        assert_eq!(g.n(), 1);
        assert_eq!(g.m(), 0);
    }

    #[test]
    fn test_path_graph() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1).unwrap();
        g.add_edge(1, 2).unwrap();
        assert_eq!(g.m(), 2);
        assert_eq!(g.degree(1), 2);
    }

    #[test]
    fn test_duplicate_edge() {
        let mut g = Graph::new(2);
        g.add_edge(0, 1).unwrap();
        g.add_edge(0, 1).unwrap();
        assert_eq!(g.m(), 1);
    }

    #[test]
    fn test_self_loop_ignored() {
        let mut g = Graph::new(2);
        g.add_edge(0, 0).unwrap();
        assert_eq!(g.m(), 0);
    }

    #[test]
    fn test_out_of_bounds() {
        let mut g = Graph::new(2);
        assert!(g.add_edge(0, 5).is_err());
    }

    #[test]
    fn test_laplacian_path3() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1).unwrap();
        g.add_edge(1, 2).unwrap();
        let l = g.laplacian();
        assert_eq!(l[0][0], 1.0);
        assert_eq!(l[1][1], 2.0);
        assert_eq!(l[0][1], -1.0);
        assert_eq!(l[2][1], -1.0);
    }

    #[test]
    fn test_connected() {
        let mut g = Graph::new(4);
        g.add_edge(0, 1).unwrap();
        g.add_edge(1, 2).unwrap();
        g.add_edge(2, 3).unwrap();
        assert!(g.is_connected());
    }

    #[test]
    fn test_disconnected() {
        let mut g = Graph::new(4);
        g.add_edge(0, 1).unwrap();
        g.add_edge(2, 3).unwrap();
        assert!(!g.is_connected());
    }

    #[test]
    fn test_fiedler_path4() {
        // Path 0-1-2-3: λ₂ ≈ 0.5858 (= 2 - √2)
        let mut g = Graph::new(4);
        g.add_edge(0, 1).unwrap();
        g.add_edge(1, 2).unwrap();
        g.add_edge(2, 3).unwrap();
        let (lambda, vec) = g.fiedler().unwrap();
        assert!((lambda - (2.0 - 2_f64.sqrt())).abs() < 0.01, "got lambda={lambda}");
        // Fiedler vector of a path has one sign change
        let pos = vec.iter().filter(|&&x| x > 0.01).count();
        let neg = vec.iter().filter(|&&x| x < -0.01).count();
        assert!(pos > 0 && neg > 0);
    }

    #[test]
    fn test_fiedler_disconnected() {
        let mut g = Graph::new(4);
        g.add_edge(0, 1).unwrap();
        g.add_edge(2, 3).unwrap();
        assert!(g.fiedler().is_err());
    }

    #[test]
    fn test_eigenmodes() {
        let mut g = Graph::new(4);
        for i in 0..3 {
            g.add_edge(i, i + 1).unwrap();
        }
        let modes = g.eigenmodes(2).unwrap();
        assert_eq!(modes.len(), 2);
        assert!(modes[0].0 < modes[1].0);
    }

    #[test]
    fn test_complete_graph_fiedler() {
        // K_n has λ₂ = n
        let n = 5;
        let mut g = Graph::new(n);
        for i in 0..n {
            for j in (i + 1)..n {
                g.add_edge(i, j).unwrap();
            }
        }
        let (lambda, _) = g.fiedler().unwrap();
        assert!((lambda - n as f64).abs() < 0.1, "got lambda={lambda}");
    }

    #[test]
    fn test_cycle_graph_fiedler() {
        // C_n has λ₂ = 2 - 2cos(2π/n)
        let n = 6;
        let mut g = Graph::new(n);
        for i in 0..n {
            g.add_edge(i, (i + 1) % n).unwrap();
        }
        let (lambda, _) = g.fiedler().unwrap();
        let expected = 2.0 - 2.0 * (2.0 * std::f64::consts::PI / n as f64).cos();
        assert!((lambda - expected).abs() < 0.05, "got {lambda}, expected {expected}");
    }

    #[test]
    fn test_star_graph_fiedler() {
        // Star K_{1,4}: λ₂ = 1
        let mut g = Graph::new(5);
        for i in 1..5 {
            g.add_edge(0, i).unwrap();
        }
        let (lambda, _) = g.fiedler().unwrap();
        assert!((lambda - 1.0).abs() < 0.05, "got {lambda}");
    }
}
