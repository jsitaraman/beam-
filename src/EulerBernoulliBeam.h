namespace beam_pp {
class EulerBernoulliBeam: public beam
{
public:
  EulerBernoulliBeam() {};
  ~EulerBernoulliBeam() {};
  EulerBernoulliBeam(std::unordered_map<std::string, std::vector<double>> & input_data);
  void Assemble(void);
private:
  const int kMaxdof_=14;     // max degrees of freedom per element
  const int nprops_max_ = 11; 
  void centrifugal_stiffening();
};
} // namespace beam_pp
