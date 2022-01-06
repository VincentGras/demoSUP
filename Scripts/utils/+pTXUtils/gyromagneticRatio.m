function g = gyromagneticRatio(nucleus)
  switch lower (nucleus)
      case {'1h', 'proton'}
          g = +42.577e6;
      otherwise
          error('unknown gyromagnetic ratio');
  end