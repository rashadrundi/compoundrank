const configuredBaseUrl =
  import.meta.env.VITE_API_BASE_URL ?? '/api/v1';

export const API_BASE_URL = configuredBaseUrl.replace(/\/$/, '');

export class ApiError extends Error {
  constructor(
    public readonly status: number,
    message: string,
    public readonly details?: unknown,
  ) {
    super(message);
    this.name = 'ApiError';
  }
}

async function readResponseBody(
  response: Response,
): Promise<unknown> {
  const text = await response.text();

  if (!text) {
    return undefined;
  }

  const contentType = response.headers.get('content-type');

  if (contentType?.includes('application/json')) {
    try {
      return JSON.parse(text);
    } catch {
      return text;
    }
  }

  return text;
}

export async function apiRequest<T>(
  path: string,
  options: RequestInit = {},
): Promise<T> {
  const normalizedPath = path.startsWith('/')
    ? path
    : `/${path}`;

  const headers = new Headers(options.headers);

  if (!headers.has('Accept')) {
    headers.set('Accept', 'application/json');
  }

  const response = await fetch(
    `${API_BASE_URL}${normalizedPath}`,
    {
      ...options,
      headers,
    },
  );

  const responseBody = await readResponseBody(response);

  if (!response.ok) {
    let message =
      `Request failed with status ${response.status}`;

    if (
      responseBody &&
      typeof responseBody === 'object' &&
      'detail' in responseBody &&
      typeof responseBody.detail === 'string'
    ) {
      message = responseBody.detail;
    }

    throw new ApiError(
      response.status,
      message,
      responseBody,
    );
  }

  return responseBody as T;
}