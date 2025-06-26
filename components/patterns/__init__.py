"""
Patterns 모듈 패키지 초기화
SMARTS 패턴 생성 및 매칭 관련 기능들을 모아놓은 패키지
"""

from .smarts_generator import SmartsGenerator
from .pattern_matcher import PatternMatcher

__all__ = [
    'SmartsGenerator',
    'PatternMatcher'
]

__version__ = '1.0.0'
__author__ = 'Chemical Structure Editor Team'
__description__ = 'SMARTS pattern generation and matching components'
