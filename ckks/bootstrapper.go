package ckks

import (
	"encoding/gob"
	"fmt"
	"math"
	"math/cmplx"
	"os"
	"time"

	// "github.com/dwkim606/test_lattigo/ckks/bettersine"
	// "github.com/dwkim606/test_lattigo/rlwe"
	// "github.com/dwkim606/test_lattigo/utils"
	"github.com/supporter-park/lattigo_optimalconv_serialize/ckks/bettersine"
	"github.com/supporter-park/lattigo_optimalconv_serialize/ring"
	"github.com/supporter-park/lattigo_optimalconv_serialize/rlwe"
	"github.com/supporter-park/lattigo_optimalconv_serialize/utils"
)

type NonPointerBtpMembers struct {
	params Parameters

	dslots    int // Number of plaintext slots after the re-encoding
	logdslots int

	encoder encoderComplex128Literal // Encoder

	prescale  float64 // Q[0]/(Q[0]/|m|)
	postscale float64 // Qi sineeval/(Q[0]/|m|)
	sinescale float64 // Qi sineeval
	sqrt2pi   float64 // (1/2pi)^{-2^r}
	scFac     float64 // 2^{r}

	coeffsToSlotsDiffScale complex128 // Matrice rescaling
	slotsToCoeffsDiffScale complex128 // Matrice rescaling

	rotKeyIndex []int // a list of the required rotation keys
}

func (btp *Bootstrapper) build_NPBM() (npbm NonPointerBtpMembers) {

	npbm.dslots = btp.dslots
	npbm.logdslots = btp.logdslots

	// npbm.encoder = GetLiteralEncoder(btp.encoder)

	npbm.prescale = btp.prescale
	npbm.postscale = btp.postscale
	npbm.sinescale = btp.sinescale
	npbm.sqrt2pi = btp.sqrt2pi
	npbm.scFac = btp.scFac

	npbm.coeffsToSlotsDiffScale = btp.coeffsToSlotsDiffScale
	npbm.slotsToCoeffsDiffScale = btp.slotsToCoeffsDiffScale

	npbm.rotKeyIndex = btp.rotKeyIndex

	return npbm
}

type PointerBtpMembers struct {
	sineEvalPoly ChebyshevInterpolation
	arcSinePoly  Poly

	pDFT    []ModdedPtDiagMatrix
	pDFTInv []ModdedPtDiagMatrix
}

type ModdedPtDiagMatrix struct {
	LogSlots   int                  // Log of the number of slots of the plaintext (needed to compute the appropriate rotation keys)
	N1         int                  // N1 is the number of inner loops of the baby-step giant-step algo used in the evaluation.
	Level      int                  // Level is the level at which the matrix is encoded (can be circuit dependant)
	Scale      float64              // Scale is the scale at which the matrix is encoded (can be circuit dependant)
	Vec        map[int][2]ring.Poly // Vec is the matrix, in diagonal form, where each entry of vec is an indexed non zero diagonal.
	naive      bool
	isGaussian bool // Each diagonal of the matrix is of the form [k, ..., k] for k a gaussian integer
}

func (btp *Bootstrapper) build_PBM() (pbm PointerBtpMembers) {

	pbm.sineEvalPoly = *btp.sineEvalPoly
	pbm.arcSinePoly = *btp.arcSinePoly

	pbm.pDFT = []ModdedPtDiagMatrix{}
	for idx := range btp.pDFT {
		var tmp ModdedPtDiagMatrix
		tmp.LogSlots = btp.pDFT[idx].LogSlots
		tmp.N1 = btp.pDFT[idx].N1
		tmp.Level = btp.pDFT[idx].Level
		tmp.Scale = btp.pDFT[idx].Scale
		tmp.naive = btp.pDFT[idx].naive
		tmp.isGaussian = btp.pDFT[idx].isGaussian

		for k, v := range btp.pDFT[idx].Vec {
			tmp.Vec[k] = [2]ring.Poly{*v[0], *v[1]}
		}

		pbm.pDFT = append(pbm.pDFT, tmp)
	}
	pbm.pDFTInv = []ModdedPtDiagMatrix{}
	for idx := range btp.pDFTInv {
		var tmp ModdedPtDiagMatrix
		tmp.LogSlots = btp.pDFTInv[idx].LogSlots
		tmp.N1 = btp.pDFTInv[idx].N1
		tmp.Level = btp.pDFTInv[idx].Level
		tmp.Scale = btp.pDFTInv[idx].Scale
		tmp.naive = btp.pDFTInv[idx].naive
		tmp.isGaussian = btp.pDFTInv[idx].isGaussian

		for k, v := range btp.pDFTInv[idx].Vec {
			tmp.Vec[k] = [2]ring.Poly{*v[0], *v[1]}
		}

		pbm.pDFTInv = append(pbm.pDFTInv, tmp)
	}

	return pbm
}

func (btp *Bootstrapper) offload_btp(npbm NonPointerBtpMembers, pbm PointerBtpMembers) {

	err := os.MkdirAll("./btp", os.ModePerm)
	if err != nil {
		fmt.Printf("Fail to make directory: %v\n", err)
	} else {
		fmt.Printf("Directory '%s' is already exist or succesfully created.\n", "./btp")
	}

	npbmFile, _ := os.Create("./btp/" + btp.identifier + "npbm.gob")
	npbmEncoder := gob.NewEncoder(npbmFile)
	err = npbmEncoder.Encode(npbm)
	if err != nil {
		fmt.Println("NPBM encoding error:", err)
		return
	}

	pbmFile, _ := os.Create("./btp/" + btp.identifier + "pbm.gob")
	pbmEncoder := gob.NewEncoder(pbmFile)
	err = pbmEncoder.Encode(pbm)
	if err != nil {
		fmt.Println("PBM encoding error:", err)
		return
	}
}

func onload_btp(npbmFileName string, pbmFileName string) (btp Bootstrapper) {

	var npbm NonPointerBtpMembers
	var pbm PointerBtpMembers

	npbmFile, _ := os.Open(npbmFileName)
	defer npbmFile.Close()
	npbmDecoder := gob.NewDecoder(npbmFile)
	err := npbmDecoder.Decode(&npbm)
	if err != nil {
		fmt.Println("NPBM decoding error:", err)
		return
	}

	pbmFile, _ := os.Open(pbmFileName)
	defer pbmFile.Close()
	pbmDecoder := gob.NewDecoder(pbmFile)
	err = pbmDecoder.Decode(&pbm)
	if err != nil {
		fmt.Println("PBM decoding error:", err)
		return
	}

	// Non-pointer things

	return
}

// Bootstrapper is a struct to stores a memory pool the plaintext matrices
// the polynomial approximation and the keys for the bootstrapping.
type Bootstrapper struct {
	*evaluator
	BootstrappingParameters
	*BootstrappingKey
	params Parameters

	dslots    int // Number of plaintext slots after the re-encoding
	logdslots int

	encoder Encoder // Encoder

	prescale     float64                 // Q[0]/(Q[0]/|m|)
	postscale    float64                 // Qi sineeval/(Q[0]/|m|)
	sinescale    float64                 // Qi sineeval
	sqrt2pi      float64                 // (1/2pi)^{-2^r}
	scFac        float64                 // 2^{r}
	sineEvalPoly *ChebyshevInterpolation // Coefficients of the Chebyshev Interpolation of sin(2*pi*x) or cos(2*pi*x/r)
	arcSinePoly  *Poly                   // Coefficients of the Taylor series of arcsine(x)

	coeffsToSlotsDiffScale complex128      // Matrice rescaling
	slotsToCoeffsDiffScale complex128      // Matrice rescaling
	pDFT                   []*PtDiagMatrix // Matrice vectors
	pDFTInv                []*PtDiagMatrix // Matrice vectors

	rotKeyIndex []int // a list of the required rotation keys

	identifier string
}

func sin2pi2pi(x complex128) complex128 {
	return cmplx.Sin(6.283185307179586*x) / 6.283185307179586
}

func cos2pi(x complex128) complex128 {
	return cmplx.Cos(6.283185307179586 * x)
}

// NewBootstrapper creates a new Bootstrapper.
func NewBootstrapper(params Parameters, btpParams *BootstrappingParameters, btpKey BootstrappingKey) (btp *Bootstrapper, err error) {

	if btpParams.SinType == SinType(Sin) && btpParams.SinRescal != 0 {
		return nil, fmt.Errorf("cannot use double angle formul for SinType = Sin -> must use SinType = Cos")
	}

	btp = newBootstrapper(params, btpParams)

	btp.BootstrappingKey = &BootstrappingKey{btpKey.Rlk, btpKey.Rtks}
	if err = btp.CheckKeys(); err != nil {
		return nil, fmt.Errorf("invalid bootstrapping key: %w", err)
	}
	btp.evaluator = btp.evaluator.WithKey(rlwe.EvaluationKey{Rlk: btpKey.Rlk, Rtks: btpKey.Rtks}).(*evaluator)

	return btp, nil
}

// NewBootstrapper creates a new Bootstrapper.
// modified by DWKIM to set coefficients multiplied on StoC matrices to be 1.0
func NewBootstrapper_mod(params Parameters, btpParams *BootstrappingParameters, btpKey BootstrappingKey) (btp *Bootstrapper, err error) {

	if btpParams.SinType == SinType(Sin) && btpParams.SinRescal != 0 {
		return nil, fmt.Errorf("cannot use double angle formul for SinType = Sin -> must use SinType = Cos")
	}

	btp = newBootstrapper(params, btpParams)
	btp.pDFT = btp.BootstrappingParameters.GenSlotsToCoeffsMatrix(1.0, btp.encoder)

	btp.BootstrappingKey = &BootstrappingKey{btpKey.Rlk, btpKey.Rtks}
	if err = btp.CheckKeys(); err != nil {
		return nil, fmt.Errorf("invalid bootstrapping key: %w", err)
	}
	btp.evaluator = btp.evaluator.WithKey(rlwe.EvaluationKey{Rlk: btpKey.Rlk, Rtks: btpKey.Rtks}).(*evaluator)

	return btp, nil
}

func NewBootstrapper_time(params Parameters, btpParams *BootstrappingParameters, btpKey BootstrappingKey) (btp *Bootstrapper, err error) {

	if btpParams.SinType == SinType(Sin) && btpParams.SinRescal != 0 {
		return nil, fmt.Errorf("cannot use double angle formul for SinType = Sin -> must use SinType = Cos")
	}

	btp = newBootstrapper_time(params, btpParams)
	btp.pDFT = btp.BootstrappingParameters.GenSlotsToCoeffsMatrix(1.0, btp.encoder)

	btp.BootstrappingKey = &BootstrappingKey{btpKey.Rlk, btpKey.Rtks}
	if err = btp.CheckKeys(); err != nil {
		return nil, fmt.Errorf("invalid bootstrapping key: %w", err)
	}
	btp.evaluator = btp.evaluator.WithKey(rlwe.EvaluationKey{Rlk: btpKey.Rlk, Rtks: btpKey.Rtks}).(*evaluator)

	return btp, nil
}

func NewBootstrapper_v8(params Parameters, btpParams *BootstrappingParameters, btpKey BootstrappingKey, id string) (btp *Bootstrapper, err error) {

	if btpParams.SinType == SinType(Sin) && btpParams.SinRescal != 0 {
		return nil, fmt.Errorf("cannot use double angle formul for SinType = Sin -> must use SinType = Cos")
	}

	btp = newBootstrapper(params, btpParams)
	btp.pDFT = btp.BootstrappingParameters.GenSlotsToCoeffsMatrix(1.0, btp.encoder)

	btp.BootstrappingKey = &BootstrappingKey{btpKey.Rlk, btpKey.Rtks}
	if err = btp.CheckKeys(); err != nil {
		return nil, fmt.Errorf("invalid bootstrapping key: %w", err)
	}
	btp.evaluator = btp.evaluator.WithKey(rlwe.EvaluationKey{Rlk: btpKey.Rlk, Rtks: btpKey.Rtks}).(*evaluator)
	btp.identifier = id

	return btp, nil
}

// newBootstrapper is a constructor of "dummy" bootstrapper to enable the generation of bootstrapping-related constants
// without providing a bootstrapping key. To be replaced by a proper factorization of the bootstrapping pre-computations.
func newBootstrapper(params Parameters, btpParams *BootstrappingParameters) (btp *Bootstrapper) {
	btp = new(Bootstrapper)

	btp.params = params
	btp.BootstrappingParameters = *btpParams.Copy()

	btp.dslots = params.Slots()
	btp.logdslots = params.LogSlots()
	if params.LogSlots() < params.MaxLogSlots() {
		btp.dslots <<= 1
		btp.logdslots++
	}

	btp.prescale = math.Exp2(math.Round(math.Log2(float64(params.Q()[0]) / btp.MessageRatio)))
	btp.sinescale = math.Exp2(math.Round(math.Log2(btp.SineEvalModuli.ScalingFactor)))
	btp.postscale = btp.sinescale / btp.MessageRatio

	btp.encoder = NewEncoder(params)
	btp.evaluator = NewEvaluator(params, rlwe.EvaluationKey{}).(*evaluator) // creates an evaluator without keys for genDFTMatrices

	btp.genSinePoly()
	btp.genDFTMatrices()

	btp.ctxpool = NewCiphertext(params, 1, params.MaxLevel(), 0)

	return btp
}

func newBootstrapper_time(params Parameters, btpParams *BootstrappingParameters) (btp *Bootstrapper) {
	fmt.Println("Bootstrapper generation with time measure")
	btp = new(Bootstrapper)

	btp.params = params
	btp.BootstrappingParameters = *btpParams.Copy()

	btp.dslots = params.Slots()
	btp.logdslots = params.LogSlots()
	if params.LogSlots() < params.MaxLogSlots() {
		btp.dslots <<= 1
		btp.logdslots++
	}

	btp.prescale = math.Exp2(math.Round(math.Log2(float64(params.Q()[0]) / btp.MessageRatio)))
	btp.sinescale = math.Exp2(math.Round(math.Log2(btp.SineEvalModuli.ScalingFactor)))
	btp.postscale = btp.sinescale / btp.MessageRatio

	start := time.Now()
	btp.encoder = NewEncoder(params)
	fmt.Println("Encoder generation (ms): ", time.Since(start))
	start = time.Now()
	btp.evaluator = NewEvaluator(params, rlwe.EvaluationKey{}).(*evaluator) // creates an evaluator without keys for genDFTMatrices
	fmt.Println("Evaluator generation (ms): ", time.Since(start))

	start = time.Now()
	btp.genSinePoly()
	fmt.Println("SinePoly generation (ms): ", time.Since(start))
	start = time.Now()
	btp.genDFTMatrices()
	fmt.Println("DFTMatrices generation (ms): ", time.Since(start))

	btp.ctxpool = NewCiphertext(params, 1, params.MaxLevel(), 0)

	return btp
}

// CheckKeys checks if all the necessary keys are present
func (btp *Bootstrapper) CheckKeys() (err error) {

	if btp.Rlk == nil {
		return fmt.Errorf("relinearization key is nil")
	}

	if btp.Rtks == nil {
		return fmt.Errorf("rotation key is nil")
	}

	rotMissing := []int{}
	for _, i := range btp.rotKeyIndex {
		galEl := btp.params.GaloisElementForColumnRotationBy(int(i))
		if _, generated := btp.Rtks.Keys[galEl]; !generated {
			rotMissing = append(rotMissing, i)
		}
	}

	if len(rotMissing) != 0 {
		return fmt.Errorf("rotation key(s) missing: %d", rotMissing)
	}

	return nil
}

// AddMatrixRotToList adds the rotations neede to evaluate pVec to the list rotations
func AddMatrixRotToList(pVec *PtDiagMatrix, rotations []int, slots int, repack bool) []int {

	if pVec.naive {
		for j := range pVec.Vec {
			if !utils.IsInSliceInt(j, rotations) {
				rotations = append(rotations, j)
			}
		}
	} else {
		var index int
		for j := range pVec.Vec {

			N1 := pVec.N1

			index = ((j / N1) * N1)

			if repack {
				// Sparse repacking, occurring during the first DFT matrix of the CoeffsToSlots.
				index &= 2*slots - 1
			} else {
				// Other cases
				index &= slots - 1
			}

			if index != 0 && !utils.IsInSliceInt(index, rotations) {
				rotations = append(rotations, index)
			}

			index = j & (N1 - 1)

			if index != 0 && !utils.IsInSliceInt(index, rotations) {
				rotations = append(rotations, index)
			}
		}
	}

	return rotations
}

func (btp *Bootstrapper) genDFTMatrices() {

	a := real(btp.sineEvalPoly.a)
	b := real(btp.sineEvalPoly.b)
	n := float64(btp.params.N())
	qDiff := float64(btp.params.Q()[0]) / math.Exp2(math.Round(math.Log2(float64(btp.params.Q()[0]))))

	// Change of variable for the evaluation of the Chebyshev polynomial + cancelling factor for the DFT and SubSum + evantual scaling factor for the double angle formula
	btp.coeffsToSlotsDiffScale = complex(math.Pow(2.0/((b-a)*n*btp.scFac*qDiff), 1.0/float64(btp.CtSDepth(false))), 0)

	// Rescaling factor to set the final ciphertext to the desired scale
	btp.slotsToCoeffsDiffScale = complex(math.Pow((qDiff*btp.params.Scale())/btp.postscale, 1.0/float64(btp.StCDepth(false))), 0)

	// CoeffsToSlots vectors
	btp.pDFTInv = btp.BootstrappingParameters.GenCoeffsToSlotsMatrix(btp.coeffsToSlotsDiffScale, btp.encoder)

	// SlotsToCoeffs vectors
	btp.pDFT = btp.BootstrappingParameters.GenSlotsToCoeffsMatrix(btp.slotsToCoeffsDiffScale, btp.encoder)

	// List of the rotation key values to needed for the bootstrapp
	btp.rotKeyIndex = []int{}

	//SubSum rotation needed X -> Y^slots rotations
	for i := btp.params.LogSlots(); i < btp.params.MaxLogSlots(); i++ {
		if !utils.IsInSliceInt(1<<i, btp.rotKeyIndex) {
			btp.rotKeyIndex = append(btp.rotKeyIndex, 1<<i)
		}
	}

	// Coeffs to Slots rotations
	for _, pVec := range btp.pDFTInv {
		btp.rotKeyIndex = AddMatrixRotToList(pVec, btp.rotKeyIndex, btp.params.Slots(), false)
	}

	// Slots to Coeffs rotations
	for i, pVec := range btp.pDFT {
		btp.rotKeyIndex = AddMatrixRotToList(pVec, btp.rotKeyIndex, btp.params.Slots(), (i == 0) && (btp.params.LogSlots() < btp.params.MaxLogSlots()))
	}
}

func (btp *Bootstrapper) genSinePoly() {

	K := int(btp.SinRange)
	deg := int(btp.SinDeg)
	btp.scFac = float64(int(1 << btp.SinRescal))

	if btp.ArcSineDeg > 0 {
		btp.sqrt2pi = 1.0

		coeffs := make([]complex128, btp.ArcSineDeg+1)

		coeffs[1] = 0.15915494309189535

		for i := 3; i < btp.ArcSineDeg+1; i += 2 {

			coeffs[i] = coeffs[i-2] * complex(float64(i*i-4*i+4)/float64(i*i-i), 0)

		}

		btp.arcSinePoly = NewPoly(coeffs)

	} else {
		btp.sqrt2pi = math.Pow(0.15915494309189535, 1.0/btp.scFac)
	}

	if btp.SinType == Sin {

		btp.sineEvalPoly = Approximate(sin2pi2pi, -complex(float64(K)/btp.scFac, 0), complex(float64(K)/btp.scFac, 0), deg)

	} else if btp.SinType == Cos1 {

		btp.sineEvalPoly = new(ChebyshevInterpolation)

		btp.sineEvalPoly.coeffs = bettersine.Approximate(K, deg, btp.MessageRatio, int(btp.SinRescal))

		btp.sineEvalPoly.maxDeg = btp.sineEvalPoly.Degree()
		btp.sineEvalPoly.a = complex(float64(-K)/btp.scFac, 0)
		btp.sineEvalPoly.b = complex(float64(K)/btp.scFac, 0)
		btp.sineEvalPoly.lead = true

	} else if btp.SinType == Cos2 {

		btp.sineEvalPoly = Approximate(cos2pi, -complex(float64(K)/btp.scFac, 0), complex(float64(K)/btp.scFac, 0), deg)

	} else {
		panic("Bootstrapper -> invalid sineType")
	}

	for i := range btp.sineEvalPoly.coeffs {
		btp.sineEvalPoly.coeffs[i] *= complex(btp.sqrt2pi, 0)
	}
}
