package ring

import (
	"github.com/dwkim606/test_lattigo/utils"
)

const precision = uint64(56)

type baseSampler struct {
	prng     utils.PRNG
	baseRing *Ring
}
