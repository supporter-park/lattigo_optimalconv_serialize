package ring

import (
	"github.com/supporter-park/lattigo_optimalconv_serialize/utils"
)

const precision = uint64(56)

type baseSampler struct {
	prng     utils.PRNG
	baseRing *Ring
}
